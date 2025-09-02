"""
evaluate the model using single source tissue on SheepGoatGTEx data
"""
import click
from train_sheepgoat_gtex import *
from src.train_utils import forward
from src.distributions import *
from src.eval_utils import *
from src.baselines import *
from tqdm import tqdm
import blitzgsea as blitz
import gseapy as gp
from sklearn.decomposition import PCA
import torch
import os
import networkx as nx
from collections import Counter


@click.command()
@click.option('--config')
@click.option('--modelfile')
@click.option('--expressfile')
@click.option('--metafile')
@click.option('--tiscolorfile')
@click.option('--mintissuenum', default=2, type=int)
@click.option('--minsamplenum', default=25, type=int)
@click.option('--traindonorfile')
@click.option('--testdonorfile')
@click.option('--valdonorfile')
@click.option('--sourcetissue')
@click.option('--resultsfile')
def main(config, modelfile, expressfile, metafile, tiscolorfile, mintissuenum, minsamplenum, traindonorfile, testdonorfile, valdonorfile, sourcetissue, resultsfile):
    # Initialise wandb
    wandb.init(config=config, mode='disabled')
    config = wandb.config
    print(config)

    # load data
    adata = GTEx_normalised_adata(expressfile=expressfile, metafile=metafile, tiscolorfile=tiscolorfile, mintissuenum=mintissuenum, minsamplenum=minsamplenum)

    collate_fn = Data.from_datalist

    # Split train/val/test
    donors = adata.obs['Participant ID'].values
    train_donors = np.loadtxt(traindonorfile, dtype=str)
    test_donors = np.loadtxt(testdonorfile, dtype=str)
    val_donors = np.loadtxt(valdonorfile, dtype=str)
    train_mask = np.isin(donors, train_donors)
    test_mask = np.isin(donors, test_donors)
    val_mask = np.isin(donors, val_donors)
    print(f'Number of training donors={len(train_donors)}, Number of validation donors={len(val_donors)}, Number of testing donors={len(test_donors)}')

    # Use certain GPU
    device = torch.device("cuda:{}".format(config.gpu) if torch.cuda.is_available() else "cpu")

    # Select dynamic/static node types
    config.update({'static_node_types': {'Tissue': (len(adata.obs['Tissue_idx'].unique()), config.d_tissue), 'metagenes': (config.meta_G, config.d_gene)}}, allow_val_change=True)
    config.update({'dynamic_node_types': {'Participant ID': (len(adata.obs['Participant ID'].unique()), config.d_patient)}}, allow_val_change=True)

    # load Model
    config.G = adata.shape[-1]
    model = HypergraphNeuralNet(config).to(device)  # .double()
    model.load_state_dict(torch.load(modelfile, map_location=torch.device('cpu')))

    # Compare to baselines
    sample_corr = True
    model.eval()
    score_fn = pearson_correlation_score
    validate = False

    source_tissues = [sourcetissue]
    target_tissues = [t for t in adata.obs['Tissue'].unique() if t not in source_tissues]

    results_df = pd.DataFrame([], columns=['individual', 'score', 'source', 'target', 'method'])

    for tt in tqdm(target_tissues):
        # Name source and target tissues
        source_name = ', '.join(source_tissues)
        target_name = tt
        # Create datasets
        aux_train_dataset = HypergraphDataset(adata[train_mask],
                                        obs_source={'Tissue': source_tissues},
                                        obs_target={'Tissue': [tt]})
        source_donor_ids = aux_train_dataset.adata_source.obs['Participant ID']
        target_donor_ids = aux_train_dataset.adata_target.obs['Participant ID']
        assert (source_donor_ids.values == target_donor_ids.values).all()
        
        aux_val_dataset = HypergraphDataset(adata[val_mask],
                                        obs_source={'Tissue': source_tissues},
                                        obs_target={'Tissue': [tt]})
        source_donor_ids = aux_val_dataset.adata_source.obs['Participant ID']
        target_donor_ids = aux_val_dataset.adata_target.obs['Participant ID']
        assert (source_donor_ids.values == target_donor_ids.values).all()
            
        aux_test_dataset = HypergraphDataset(adata[test_mask],
                                        obs_source={'Tissue': source_tissues},
                                        obs_target={'Tissue': [tt]})
        source_donor_ids = aux_test_dataset.adata_source.obs['Participant ID']
        target_donor_ids = aux_test_dataset.adata_target.obs['Participant ID']
        assert (source_donor_ids.values == target_donor_ids.values).all()
        if len(target_donor_ids)<1:
            print(f'{tt}: Number of test donors={len(target_donor_ids)}. Insufficient donors. Skip.')
            continue
        
        # Prepare source expression data
        x_train_ = aux_train_dataset.adata_source.layers['x'].toarray()
        x_train_covs = aux_train_dataset.adata_source.obsm['Participant ID_feat'].toarray()
        x_val_ = aux_val_dataset.adata_source.layers['x'].toarray()
        x_val_covs = aux_val_dataset.adata_source.obsm['Participant ID_feat'].toarray()
        x_test_ = aux_test_dataset.adata_source.layers['x'].toarray()
        x_test_covs = aux_test_dataset.adata_source.obsm['Participant ID_feat'].toarray()
        y_train = aux_train_dataset.adata_target.layers['x'].toarray()
        y_val = aux_val_dataset.adata_target.layers['x'].toarray()
        y_test = aux_test_dataset.adata_target.layers['x'].toarray()
        # Append donor metadata
        x_train_aux = aux_train_dataset.adata_source.obsm['Participant ID_feat'].toarray()
        x_val_aux = aux_val_dataset.adata_source.obsm['Participant ID_feat'].toarray()
        x_test_aux = aux_test_dataset.adata_source.obsm['Participant ID_feat'].toarray()
        x_train = np.concatenate((x_train_, x_train_aux), axis=-1)
        x_val = np.concatenate((x_val_, x_val_aux), axis=-1)
        x_test = np.concatenate((x_test_, x_test_aux), axis=-1)
        if validate:
            x_test = x_val
            y_test = y_val
            x_test_ = x_val_
            x_test_covs = x_val_covs
        
        # Blood surrogate baseline
        sample_scores = score_fn(y_test, x_test_, sample_corr=sample_corr)
        
        # Append results
        scores = sample_scores
        df_ = pd.DataFrame({'individual': target_donor_ids.values,
                            'score': scores,
                            'source': [source_name] * len(scores),
                            'target': [target_name] * len(scores),
                            'method': [f'{source_name} surrogate'] * len(scores)})
        results_df = pd.concat([results_df, df_])
        
        # Mean baseline
        if y_train.shape[0] and sample_corr:  # Not defined when sampling units are genes
            means = y_train.mean(axis=0)
            y_test_pred = np.repeat(means[None, :], y_test.shape[0], axis=0)
            sample_scores = score_fn(y_test, y_test_pred, sample_corr=sample_corr)

            # Append results
            scores = sample_scores
            df_ = pd.DataFrame({'individual': target_donor_ids.values,
                                'score': scores,
                                'source': [source_name] * len(scores),
                                'target': [target_name] * len(scores),
                                'method': ['mean'] * len(scores)})
            results_df = pd.concat([results_df, df_])
        
            # KNN baseline
            x_train_knn = np.concatenate((x_train_, y_train), axis=-1)
            test_nans = np.full((x_test_.shape[0], y_train.shape[1]), np.nan)
            x_test_knn = np.concatenate((x_test_, test_nans), axis=-1)
            x_knn = np.concatenate((x_train_knn, x_test_knn), axis=0)

            x_knn_covs = np.concatenate((x_train_covs, x_test_covs), axis=0)
            knn_imp = impute_knn(x_knn, covariates=x_knn_covs, k=20)
            knn_imp_ = knn_imp[x_train_.shape[0]:, x_train_.shape[1]:]
            sample_scores = score_fn(y_test, knn_imp_, sample_corr=sample_corr)
        
            # Append results
            scores = sample_scores
            df_ = pd.DataFrame({'individual': target_donor_ids.values,
                                'score': scores,
                                'source': [source_name] * len(scores),
                                'target': [target_name] * len(scores),
                                'method': ['kNN'] * len(scores)})
            results_df = pd.concat([results_df, df_])
        
            # TEEBoT baseline.
            y_test_pred = PCA_linear_regression_baseline(x_train, y_train, x_test)        
            sample_scores = score_fn(y_test, y_test_pred, sample_corr=sample_corr)
        
            # Append results
            scores = sample_scores
            df_ = pd.DataFrame({'individual': target_donor_ids.values,
                                'score': scores,
                                'source': [source_name] * len(scores),
                                'target': [target_name] * len(scores),
                                'method': [f'PCA'] * len(scores)})
            results_df = pd.concat([results_df, df_])
        
            # Hypergraph baseline
            aux_train_loader = DataLoader(aux_train_dataset, batch_size=config.batch_size, collate_fn=collate_fn, shuffle=True, drop_last=True)
            aux_test_loader = DataLoader(aux_test_dataset, batch_size=len(aux_test_dataset), collate_fn=collate_fn, shuffle=False)

            # Compute predictions and score
            model.eval()
            with torch.no_grad():
                if validate:
                    d = next(iter(aux_val_loader))
                else:
                    d = next(iter(aux_test_loader))

                out, node_features = forward(d, model, device, preprocess_fn=None)
                y_test_pred = out['px_rate'].cpu().numpy()  # torch.distributions.normal.Normal(loc=out['px_rate'], scale=out['px_r']).mean.cpu().numpy()
                y_test_ = d.x_target.cpu().numpy()
            assert np.allclose(y_test_, y_test)

            sample_scores = score_fn(y_test, y_test_pred, sample_corr=sample_corr)

            # Append results
            scores = sample_scores
            df_ = pd.DataFrame({'individual': target_donor_ids.values,
                                'score': scores,
                                'source': [source_name] * len(scores),
                                'target': [target_name] * len(scores),
                                'method': [f'HYFA ({source_name})'] * len(scores)})
            results_df = pd.concat([results_df, df_])

            # Hypergraph baseline (all tissues)
            aux_test_dataset = HypergraphDataset(adata[test_mask],
                                            obs_source={'Participant ID': list(aux_test_dataset.donor_map.values()),
                                                        'Tissue': [t for t in adata.uns['Tissue_dict'].keys() if t != tt]},
                                            obs_target={'Tissue': [tt]})

            aux_test_loader = DataLoader(aux_test_dataset, batch_size=len(aux_test_dataset), collate_fn=collate_fn, shuffle=False)

            # Compute predictions and score
            model.eval()
            with torch.no_grad():
                if validate:
                    d = next(iter(aux_val_loader))
                else:
                    d = next(iter(aux_test_loader))

                out, node_features = forward(d, model, device, preprocess_fn=None)
                y_test_pred = out['px_rate'].cpu().numpy()  # torch.distributions.normal.Normal(loc=out['px_rate'], scale=out['px_r']).mean.cpu().numpy()
                y_test_ = d.x_target.cpu().numpy()

            sample_scores = score_fn(y_test, y_test_pred, sample_corr=sample_corr)

            # Append results
            scores = sample_scores
            df_ = pd.DataFrame({'individual': target_donor_ids.values,
                                'score': scores,
                                'source': [source_name] * len(scores),
                                'target': [target_name] * len(scores),
                                'method': ['HYFA (all)'] * len(scores)})
            results_df = pd.concat([results_df, df_])
        results_df.to_csv(resultsfile, index=None)


if __name__ == '__main__':
    main()
