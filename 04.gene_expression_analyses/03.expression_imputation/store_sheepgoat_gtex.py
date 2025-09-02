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


def create_dataframe(participant_ids, tissue_ids, expression, donor_map, tissue_dict_inv, symbols):
    participant_ids = np.concatenate(participant_ids, axis=0)
    tissue_ids = np.concatenate(tissue_ids, axis=0)
    expression = np.concatenate(expression, axis=0)
    df_metadata = pd.DataFrame({'Participant ID': [donor_map[p] for p in participant_ids],
                               'Tissue': [tissue_dict_inv[t] for t in tissue_ids]})
    df = pd.DataFrame(expression, columns=symbols)
    df = pd.concat([df_metadata, df], axis=1)
    df = df.set_index('Participant ID')
    return df

@click.command()
@click.option('--config')
@click.option('--modelfile')
@click.option('--expressfile')
@click.option('--metafile')
@click.option('--tiscolorfile')
@click.option('--mintissuenum', default=2, type=int)
@click.option('--minsamplenum', default=25, type=int)
@click.option('--observedfile')
@click.option('--imputedfile')
def main(config, modelfile, expressfile, metafile, tiscolorfile, mintissuenum, minsamplenum, observedfile, imputedfile):
    # Initialise wandb
    wandb.init(config=config, mode='disabled')
    config = wandb.config
    print(config)

    # load data
    adata = GTEx_normalised_adata(expressfile=expressfile, metafile=metafile, tiscolorfile=tiscolorfile, mintissuenum=mintissuenum, minsamplenum=minsamplenum)
    dataset = HypergraphDataset(adata, static=True)
    tissue_dict = adata.uns['Tissue_dict']
    tissue_dict_inv = {v: k for k, v in tissue_dict.items()}
    
    # Use GPU/CPU
    device = torch.device("cuda:{}".format(config.gpu) if torch.cuda.is_available() else "cpu")
    print(torch.cuda.is_available())
    print(device)
    
    # Select dynamic/static node types
    config.static_node_types = {'Tissue': (len(adata.obs['Tissue_idx'].unique()), config.d_tissue),
                                'metagenes': (config.meta_G, config.d_gene)}
    config.dynamic_node_types = {'Participant ID': (len(adata.obs['Participant ID'].unique()), config.d_patient)}

    # Load Model
    config.G = adata.shape[-1]
    model = HypergraphNeuralNet(config).to(device)  # .double()
    model.load_state_dict(torch.load(modelfile, map_location=torch.device('cpu')))
    model.eval()

    source_participant_ids = []
    source_tissue_ids = []
    source_expression = []
    target_participant_ids = []
    target_tissue_ids = []
    target_expression = []

    for i, d in tqdm(enumerate(dataset)):
        # Set target tissues to missing tissues
        d.target['Tissue'] = torch.tensor([t for t in np.arange(len(tissue_dict)) if t not in d.source['Tissue']])
        d.target['Participant ID'] = torch.zeros_like(d.target['Tissue']) + d.source['Participant ID'][0]
        
        # Make predictions
        with torch.no_grad():
            out, node_features = forward(d, model, device, preprocess_fn=None) 
            y_pred = out['px_rate']
        
        # Store
        source_participant_ids.append(d.source['Participant ID'].cpu().numpy() + i)
        source_tissue_ids.append(d.source['Tissue'].cpu().numpy())
        source_expression.append(d.x_source.cpu().numpy())
        target_participant_ids.append(d.target['Participant ID'].cpu().numpy() + i)
        target_tissue_ids.append(d.target['Tissue'].cpu().numpy())
        target_expression.append(y_pred.cpu().numpy())

    # Store data in dataframes
    df_imputed = create_dataframe(target_participant_ids, target_tissue_ids, target_expression,
                                donor_map=dataset.donor_map,
                                tissue_dict_inv=tissue_dict_inv,
                                symbols=adata.var['Symbol'])
    df_observed = create_dataframe(source_participant_ids, source_tissue_ids, source_expression,
                                donor_map=dataset.donor_map,
                                tissue_dict_inv=tissue_dict_inv,
                                symbols=adata.var['Symbol'])

    df_observed.to_csv(observedfile)
    df_imputed.to_csv(imputedfile)

if __name__ == '__main__':
    main()
