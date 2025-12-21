import sctour as sct
import scanpy as sc
import numpy as np
import pandas as pd
import torch
import anndata
import argparse
import gc
import matplotlib.pyplot as plt
from pathlib import Path
DATA_DIR = Path("/root/autodl-tmp/dataset")

global_parser = argparse.ArgumentParser(description="compute velocity")
global_parser.add_argument('--dataset', type=str, help="the name of the dataset")
global_parser.add_argument('--clusterkey', type=str, help="the key of cluster")
global_parser.add_argument('--X_embed', type=str, help="the key of X_embed")
global_parser.add_argument('--K_FOLD', type=int, required=True, help="the number of K_FOLD")
global_parser.add_argument('--seed', type=int, help="random seed")
global_args = global_parser.parse_args()

DATASET = global_args.dataset
CLUSTER_KEY = global_args.clusterkey
VIS_EMB = global_args.X_embed
K_FOLD = global_args.K_FOLD
SEED = global_args.seed
SAVE_DATA = True
if SAVE_DATA:
    (DATA_DIR / DATASET / "processed").mkdir(parents=True, exist_ok=True)

def run_sctour_mse(fold):
    adata = sc.read_h5ad(DATA_DIR / DATASET / "processed" / f"adata_preprocessed_{fold}.h5ad")
    adata.X = adata.layers["spliced"]
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    tnode = sct.train.Trainer(adata, 
                          loss_mode='mse', 
                          alpha_recon_lec=0.5, 
                          alpha_recon_lode=0.5,
                          random_state=SEED,
                         batch_norm=False)
    tnode.device = torch.device('cuda:0')
    tnode.train()
    adata.obs['ptime'] = tnode.get_time()
    mix_zs, zs, pred_zs = tnode.get_latentsp(alpha_z=0.5, alpha_predz=0.5)
    adata.obsm['X_TNODE'] = mix_zs
    adata.obsm['X_VF'] = tnode.get_vector_field(adata.obs['ptime'].values, adata.obsm['X_TNODE'])

    new_adata = anndata.AnnData(adata.obsm['X_TNODE'].copy())
    new_adata.layers['spliced'] = adata.obsm['X_TNODE'].copy()
    new_adata.layers['unspliced'] = adata.obsm['X_TNODE'].copy()
    new_adata.layers['sctour_mse_velocity'] = adata.obsm['X_VF'].copy()
    new_adata.obs.index = adata.obs.index.copy()
    for key in adata.obs:
        new_adata.obs[key] = adata.obs[key].copy()
    new_adata.obs['sctour_mse_time'] = new_adata.obs['ptime']
    for key in adata.uns:
        new_adata.uns[key] = adata.uns[key].copy()
    for key in adata.obsm:
        new_adata.obsm[key] = adata.obsm[key].copy()
    for key in adata.obsp:
        new_adata.obsp[key] = adata.obsp[key].copy()

    new_adata.write_h5ad(DATA_DIR / DATASET / "processed" / f"adata_run_sctour_mse_{fold}.h5ad")
    del adata, new_adata
    gc.collect()

if K_FOLD == 0:
    fold = 'full'
    run_sctour_mse(fold)
else:
    for fold in range(K_FOLD):
        run_sctour_mse(fold)
