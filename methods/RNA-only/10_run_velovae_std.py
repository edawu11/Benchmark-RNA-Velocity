import argparse
import anndata
import numpy as np
import scvelo as scv
import scanpy as sc
import sys
import torch
import os.path
import velovae as vv
import pickle as pickle
import matplotlib.pyplot as plt
import time
import pandas as pd
from os.path import exists
from pathlib import Path
import gc
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

torch.manual_seed(SEED)
torch.cuda.manual_seed(SEED)
np.random.seed(SEED)

def run_velovae_std(fold):
    adata = sc.read_h5ad(DATA_DIR / DATASET / "processed" / f"adata_preprocessed_{fold}.h5ad")
    vae = vv.VAE(adata, 
             tmax=20, 
             dim_z=5,
             device='cuda:0')
    
    vae.train(adata, plot=False, embed=VIS_EMB.split('_')[1])
    vae.save_anndata(adata, 'velovae_std', DATA_DIR / DATASET / 'processed', file_name=f"adata_run_velovae_std_{fold}.h5ad")
    if torch.cuda.is_available():
        torch.cuda.empty_cache()  # Clear PyTorch's GPU memory cache
        torch.cuda.synchronize()
    del adata  # Remove reference to AnnData object
    gc.collect()

# run not-full model
if K_FOLD == 0:
    fold = 'full'
    run_velovae_std(fold)
else:
    for fold in range(K_FOLD):
        run_velovae_std(fold)