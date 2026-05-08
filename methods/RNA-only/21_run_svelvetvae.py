import argparse
import velvetvae as vt

# general packages
import numpy as np
import pandas as pd
import torch
from scipy.sparse import issparse

# velocity packages
import scanpy as sc
import scvelo as scv
import anndata as ann

# plotting packages
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm, trange
from IPython.display import clear_output

# color palette object
# from colors import colorpalette as colpal
import os

import gc

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

def run_svelvetvae(fold):
    adata = sc.read_h5ad(DATA_DIR / DATASET / "processed" / f"adata_preprocessed_{fold}.h5ad")
    adata.layers['spliced'] = (adata.layers['spliced'].A if 
                           issparse(adata.layers['spliced']) else 
                           adata.layers['spliced'])
    adata.layers['unspliced'] = (adata.layers['unspliced'].A if 
                           issparse(adata.layers['unspliced']) else 
                           adata.layers['unspliced'])
    adata.layers['total'] = adata.layers['spliced'] + adata.layers['unspliced']
    adata.layers['spliced'] = adata.layers['spliced'].astype(np.float32)
    adata.layers['unspliced'] = adata.layers['unspliced'].astype(np.float32)
    adata.layers['total'] = adata.layers['total'].astype(np.float32)
    adata.X = adata.X.astype(np.float32)
    vt.pp.neighborhood(adata, n_neighbors=30)
    vt.ut.set_seed(SEED)      
    vt.md.Svelvet.setup_anndata(adata, x_layer='total', u_layer='unspliced', knn_layer='knn_index')
    model = vt.md.Svelvet(
        adata,
        n_latent = 50,
        linear_decoder = True,
        neighborhood_space="latent_space",
        gamma_mode = "learned",
    )
    model.setup_model(gamma_kwargs={'gamma_min':0.1,'gamma_max':1})
    model.train(
        batch_size = adata.shape[0],
        max_epochs = 1000, 
        freeze_vae_after_epochs = 200,
        constrain_vf_after_epochs = 200,
        lr=0.001,
    )
    V = model.predict_velocity()
    V = V.A if issparse(V) else V
    V = np.nan_to_num(V, nan=0, neginf=0, posinf=0)
    adata.layers['svelvetvae_velocity'] = V
    # scv.tl.velocity_graph(adata, vkey='svelvetvae_velocity')
    # scv.tl.velocity_pseudotime(adata, vkey='svelvetvae_velocity')
    # adata.obs['time_info'] = adata.obs['svelvetvae_velocity_pseudotime']
    if SAVE_DATA:
        adata.write_h5ad(DATA_DIR / DATASET / "processed" / f"adata_run_svelvetvae_{fold}.h5ad")
    if torch.cuda.is_available():
        torch.cuda.empty_cache()  # Clear PyTorch's GPU memory cache
        torch.cuda.synchronize()
    del adata  # Remove reference to AnnData object
    gc.collect()

if K_FOLD == 0:
    fold = 'full'
    run_svelvetvae(fold)
else:
    for fold in range(K_FOLD):
        run_svelvetvae(fold)    
    
    