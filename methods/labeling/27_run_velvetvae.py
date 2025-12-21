import argparse
import RNA_velocity_benchmark_github.methods.metablic_labeling_enhanced.velvetvae as vt

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

def run_velvetvae(fold):
    adata = sc.read_h5ad(DATA_DIR / DATASET / "processed" / f"adata_preprocessed_{fold}.h5ad")

    adata.layers['total'] = (
        adata.layers['total'].A if issparse(adata.layers['total']) else adata.layers['total']
    )
    adata.layers['new'] = (
        adata.layers['new'].A if issparse(adata.layers['new']) else adata.layers['new']
    )
    
    adata = vt.pp.size_normalize(
        adata, 
        genes=list(adata.var_names), 
        total_layer='total', 
        new_layer='new',
        unsparsify=True
    )
    adata.layers["total"] = adata.layers["total"].astype(np.float32)
    adata.layers["new"] = adata.layers["new"].astype(np.float32)
    
    vt.pp.neighborhood(adata, n_neighbors=30)
    
    vt.ut.set_seed(SEED)
    
    vt.md.Velvet.setup_anndata(adata, x_layer='total', n_layer='new', knn_layer='knn_index')

    model = vt.md.Velvet(
        adata,
        n_latent = 50,
        linear_decoder = True,
        neighborhood_space="latent_space",
        biophysical_model = "full",
        gamma_mode = "learned",
        labelling_time = 2.0,
    )

    model.setup_model()
    
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
    adata.layers['velvetvae_velocity'] = V
    if SAVE_DATA:
        adata.write_h5ad(DATA_DIR / DATASET / "processed" / f"adata_run_velvetvae_{fold}.h5ad")
    if torch.cuda.is_available():
        torch.cuda.empty_cache()  # Clear PyTorch's GPU memory cache
        torch.cuda.synchronize()
    del adata  # Remove reference to AnnData object
    gc.collect()

if K_FOLD == 0:
    fold = 'full'
    run_velvetvae(fold)
else:
    for fold in range(K_FOLD):
        run_velvetvae(fold)    
    
    