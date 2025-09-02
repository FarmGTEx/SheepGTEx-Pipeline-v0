"""
Trains the model on SheepGoatGTEx data
"""
import torch
import numpy as np
import pandas as pd
import os
import wandb
import argparse
from torch.utils.data import Dataset, DataLoader
from src.hnn import HypergraphNeuralNet
from src.data import Data
from src.dataset import HypergraphDataset
from src.data_utils import *
from src.eval_utils import *
from src.train_utils import train
import scanpy as sc


def GTEx(file, tissueselect, mintissuenum=0, minsamplenum=0):
    """
    Loads processed GTEx data
    :param file: path of the CSV file
    :return: Returns
        - data: numpy array of shape=(nb_samples, nb_genes)
        - gene_symbols: numpy array with gene symbols. Shape=(nb_genes,)
        - sampl_ids: numpy array with sample IDs (GTEx IDs of individuals, e.g. GTEX-1117F). Shape=(nb_samples,)
        - tissues: numpy array indicating the tissue of each sample. Shape=(nb_samples,)
    """
    # Load data
    df = pd.read_csv(file, index_col=0)  # .sample(frac=1, random_state=random_seed)
    df.index.name = 'SUBJID'
    # select tissues
    df = df[df['tissue'].isin(tissueselect)]
    df = df.sort_values(by=['SUBJID', 'tissue'])
    # Delete under-represented participants and tissues
    if mintissuenum and minsamplenum:
        print(f'mintissuenum={mintissuenum}, minsamplenum={minsamplenum}')
        for i in range(5):
            df = df.loc[df.index.value_counts()>=mintissuenum]
            tiscounts = df['tissue'].value_counts()
            df = df[df['tissue'].isin(tiscounts[tiscounts>=minsamplenum].index)]
        print(f'Number of samples={df.shape[0]}, Number of donors={len(df.index.unique())}, Number of tissues={len(df.tissue.unique())}')
    tissues = df['tissue'].values
    sampl_ids = df.index.values
    del df['tissue']
    data = np.float32(df.values)
    gene_symbols = df.columns.values
    return data, gene_symbols, sampl_ids, tissues

def GTEx_metadata(file):
    """
    Loads metadata DataFrame with information about individuals
    :param file: path of the file
    :return: Pandas DataFrame with subjects' information
    """
    df = pd.read_csv(file, delimiter='\t')
    df = df.set_index('SUBJID')
    df = df.sort_index()
    return df

def GTEx_normalised_adata(expressfile, metafile, tiscolorfile, mintissuenum=0, minsamplenum=0):
    # Set up tissue colors
    with open(tiscolorfile,'r') as f:
        tissue2color = {}
        for line in f:
            key, value = line.strip().split()
            tissue2color[key] = value

    data, gene_symbols, sampl_ids, tissues = GTEx(file=expressfile, tissueselect=tissue2color.keys(), mintissuenum=mintissuenum, minsamplenum=minsamplenum)
    metadata_df = GTEx_metadata(file=metafile)

    adata = sc.AnnData(data)
    adata.var['Symbol'] = gene_symbols
    adata.obs['Participant ID'] = sampl_ids
    adata.obs['Tissue'] = tissues

    # Static keys
    adata.obs['Tissue_idx'], tissue_dict = map_to_ids(adata.obs['Tissue'].values)
    adata.uns['Tissue_dict'] = tissue_dict
    # del adata.obs['Tissue']

    # Dynamic keys
    adata.obs['Participant ID_dyn'] = adata.obs['Participant ID']

    # Populate participant features
    adata.obs['Age'] = [float(a) for a in metadata_df.loc[adata.obs['Participant ID']]['AGE'].values]
    adata.obs['Sex'] = metadata_df.loc[adata.obs['Participant ID']]['SEX'].values
    donor_age = adata.obs['Age']
    donor_sex, donor_sex_dict = map_to_ids(adata.obs['Sex'])
    adata.obsm['Participant ID_feat'] = np.stack((donor_age, donor_sex), axis=-1)
    adata.uns['Sex_dict'] = donor_sex_dict

    # Put gene expression in layer
    adata.layers['x'] = adata.X

    tissue_dict_inv = {v: k for k, v in adata.uns['Tissue_dict'].items()}
    adata.uns['Tissue_colors'] = []
    for i in range(len(tissue_dict_inv)):
        adata.uns['Tissue_colors'].append(tissue2color[tissue_dict_inv[i]])

    return adata

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', dest='config', default='configs/default.yaml', type=str)
    parser.add_argument('--expressfile', dest='expressfile', default='sheep/express.csv', type=str)
    parser.add_argument('--metafile', dest='metafile', default='sheep/meta.txt', type=str)
    parser.add_argument('--tiscolorfile', dest='tiscolorfile', default='sheep/tissue_color.list', type=str,  help="file of selected tissues and the color code, seperated by tab")
    parser.add_argument('--mintissuenum', dest='mintissuenum', default=2, type=int, help="minimum number of tissues in a individual")
    parser.add_argument('--minsamplenum', dest='minsamplenum', default=25, type=int, help="minimum number of individuals in a tissue")
    parser.add_argument('--traindonorfile', dest='traindonorfile', default="sheep/splits/train.txt", type=str)
    parser.add_argument('--valdonorfile', dest='valdonorfile', default="sheep/splits/val.txt", type=str)
    parser.add_argument('--testdonorfile', dest='testdonorfile', default="sheep/splits/test.txt", type=str)
    parser.add_argument('--modelfile', dest='modelfile', default='sheep/sheep.model.pth', type=str)
    args, unknown = parser.parse_known_args()

    np.random.seed(0)
    num_workers = 4
    
    # Initialise wandb
    os.environ["WANDB__SERVICE_WAIT"] = "300"
    run = wandb.init(config=args.config, mode='offline')
    config = wandb.config
    print(config)

    # Load data
    adata = GTEx_normalised_adata(expressfile=args.expressfile, metafile=args.metafile, tiscolorfile=args.tiscolorfile, mintissuenum=args.mintissuenum, minsamplenum=args.minsamplenum)

    # Split train/val/test
    donors = adata.obs['Participant ID'].values
    train_donors, test_donors = split_patient_train_test(donors, train_rate=0.8)
    train_donors, val_donors = split_patient_train_test(train_donors, train_rate=0.75)
    train_mask = np.isin(donors, train_donors)
    val_mask = np.isin(donors, val_donors)
    test_mask = np.isin(donors, test_donors)
    print(f'Number of training donors={len(train_donors)}, Number of validation donors={len(val_donors)}, Number of testing donors={len(test_donors)}')

    np.savetxt(args.traindonorfile, train_donors, fmt='%s')
    np.savetxt(args.valdonorfile, val_donors, fmt='%s')
    np.savetxt(args.testdonorfile, test_donors, fmt='%s')

    collate_fn = Data.from_datalist
    dtype = torch.float32  # torch.double

    train_dataset = HypergraphDataset(adata[train_mask], dtype=dtype, disjoint=True, static=False)
    val_dataset = HypergraphDataset(adata[val_mask], dtype=dtype, disjoint=False, static=True)
    # test_dataset = HypergraphDataset(adata[test_mask], dtype=dtype, static=True)

    train_loader = DataLoader(train_dataset, batch_size=config.batch_size, collate_fn=collate_fn, shuffle=True, num_workers=num_workers)
    val_loader = DataLoader(val_dataset, batch_size=config.batch_size, collate_fn=collate_fn, shuffle=False, num_workers=num_workers)
    # test_loader = DataLoader(test_dataset, batch_size=config.batch_size, collate_fn=collate_fn, shuffle=False, num_workers=num_workers)

    # Use GPU/CPU
    device = torch.device("cuda:{}".format(config.gpu) if torch.cuda.is_available() else "cpu")
    print(torch.cuda.is_available())
    print(device)
    
    # Select dynamic/static node types
    config.static_node_types = {'Tissue': (len(adata.obs['Tissue_idx'].unique()), config.d_tissue),
                                'metagenes': (config.meta_G, config.d_gene)}
    config.dynamic_node_types = {'Participant ID': (len(adata.obs['Participant ID'].unique()), config.d_patient)}

    # Model
    config.G = adata.shape[-1]
    model = HypergraphNeuralNet(config).to(device)  # .double()

    # Train
    def rho(x, out):
        x_pred = out['px_rate'].detach().cpu().numpy()
        return np.mean(pearson_correlation_score(x, x_pred, sample_corr=True))
    metric_fns = [rho]
    train(config,
        model=model,
        loader=train_loader,
        val_loader=val_loader,
        device=device,
        preprocess_fn=None,
        compute_metrics_train=False,
        metric_fns=metric_fns)

    torch.save(model.state_dict(), args.modelfile)
    run.finish()
    
