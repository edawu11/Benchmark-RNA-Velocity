import anndata
import numpy as np
import scvelo as scv
import scanpy as sc
import sys
import torch
import os.path
import pickle as pickle
import matplotlib.pyplot as plt
import pandas as pd
import latentvelo as ltv
from os.path import exists
import gc
from pathlib import Path
import argparse
DATA_DIR = Path("/root/autodl-tmp/dataset")
MODEL_DIR = Path("/root/autodl-tmp/notebook")

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

torch.manual_seed(SEED)
torch.cuda.manual_seed(SEED)
np.random.seed(SEED)

def run_latentvelo_std(fold):
    adata = sc.read_h5ad(DATA_DIR / DATASET / "processed" / f"adata_preprocessed_{i}.h5ad")

    adata = ltv.utils.standard_clean_recipe(adata, n_top_genes = None)
    model = ltv.models.VAE(observed=adata.n_vars, latent_dim=20, zr_dim=1,h_dim=2)
    epochs, val_ae, val_traj = ltv.train(model, adata, batch_size = 128, learning_rate=1e-2,
                                         epochs=50, name=f'{DATASET}_parameters', grad_clip=100)
    gc.collect()

    latent_adata, adata = ltv.output_results(model, adata,
                                         gene_velocity = True,
										 decoded = True,
										 embedding= VIS_EMB.split('_')[1])
    latent_adata.layers["latentvelo_std_velocity"] = latent_adata.layers.pop("spliced_velocity")
    latent_adata.obs["latentvelo_std_time"] = latent_adata.obs["latent_time"]
    if SAVE_DATA:
        latent_adata.write_h5ad(DATA_DIR / DATASET / "processed" / f"adata_run_latentvelo_std_{i}.h5ad")
    if torch.cuda.is_available():
        torch.cuda.empty_cache()  # Clear PyTorch's GPU memory cache
        torch.cuda.synchronize()
    del adata  # Remove reference to AnnData object
    gc.collect()

if K_FOLD == 0:
    fold = 'full'
    run_latentvelo_std(fold)
else:
    for fold in range(K_FOLD):
        run_latentvelo_std(fold)
