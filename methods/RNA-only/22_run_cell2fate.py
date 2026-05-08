import argparse
import os
os.chdir('..')
import scvelo as scv
import scanpy as sc
import cell2fate as c2f
import pickle as pickle
from datetime import datetime
import pandas as pd
import numpy as np
from os.path import exists
import matplotlib.pyplot as plt
import torch
import time
import scipy
import torch
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

def run_cell2fate(fold):
    adata = sc.read(DATA_DIR / DATASET / "processed" / f"adata_preprocessed_{fold}.h5ad")
    if scipy.sparse.issparse(adata.layers['raw_spliced']):
        adata.layers['raw_spliced'] = np.array(adata.layers['raw_spliced'].toarray(), dtype=np.float32)
    if scipy.sparse.issparse(adata.layers['raw_unspliced']):
        adata.layers['raw_unspliced'] = np.array(adata.layers['raw_unspliced'].toarray(), dtype=np.float32)
    c2f.Cell2fate_DynamicalModel.setup_anndata(adata, spliced_label='raw_spliced', unspliced_label='raw_unspliced')
    n_modules = int(len(np.unique(adata.obs[CLUSTER_KEY]))*1.15)
    mod = c2f.Cell2fate_DynamicalModel(adata, n_modules = n_modules,
                                   Tmax_prior={"mean": 50., "sd": 50.})
    mod.train()
    adata = mod.export_posterior(adata,sample_kwargs={"num_samples": 30, "batch_size" : None,
                                                     "use_gpu" : True, 'return_samples': False})
    mod.compute_and_plot_total_velocity(adata, save = False, delete = False, plot=False)
    adata.layers["cell2fate_velocity"] = adata.layers.pop("Velocity")
    adata.layers['cell2fate_velocity'] = adata.layers['cell2fate_velocity'].numpy()
    adata.obs['cell2fate_time'] = adata.obs['Time (hours)']
    if SAVE_DATA:
        adata.write_h5ad(DATA_DIR / DATASET / "processed" / f"adata_run_cell2fate_{fold}.h5ad")
    if torch.cuda.is_available():
        torch.cuda.empty_cache()  # Clear PyTorch's GPU memory cache
        torch.cuda.synchronize()
    del adata  # Remove reference to AnnData object
    gc.collect()

if K_FOLD == 0:
    fold = 'full'
    run_cell2fate(fold)
else:
    for fold in range(K_FOLD):
        run_cell2fate(fold)













