import argparse
import anndata
import numpy as np
import scvelo as scv
import scanpy as sc
import sys
import torch
import os.path
import celldancer as cd
import pickle as pickle
import matplotlib.pyplot as plt
import pandas as pd
from os.path import exists
import celldancer.cdplt as cdplt
from celldancer.cdplt import colormap
import time
from celldancer.utilities import export_velocity_to_dynamo
import celldancer.utilities as cdutil
import gc
from pathlib import Path
import os
import shutil
METHOD = 'celldancer'
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
    (MODEL_DIR / DATASET / METHOD ).mkdir(parents=True, exist_ok=True)

def run_celldancer(fold):
    adata = sc.read_h5ad(DATA_DIR / DATASET / "processed" / f"adata_preprocessed_{fold}.h5ad")
    input_data = cdutil.adata_to_df_with_embed(adata,
                                  us_para=['Mu','Ms'],
                                  cell_type_para=CLUSTER_KEY,
                                  embed_para=VIS_EMB,
                                  save_path= MODEL_DIR / DATASET / METHOD / f'cell_type_u_s_{fold}.csv')
    loss_df, cellDancer_df=cd.velocity(input_data,
                                    gene_list= np.array(input_data['gene_name']),
                                    permutation_ratio=0.125,
                                    n_jobs=10)
    cellDancer_df=cd.compute_cell_velocity(cellDancer_df=cellDancer_df, projection_neighbor_choice='gene', 
                                        expression_scale='power10', projection_neighbor_size=10, speed_up=(100,100))
    adata = export_velocity_to_dynamo(cellDancer_df,adata)
    adata.layers["celldancer_velocity"] = adata.layers.pop("velocity_S")
    if SAVE_DATA:
        cellDancer_df.to_csv(DATA_DIR / DATASET / "processed" / f"cellDancer_df_{fold}.csv", index=False)
        adata.write_h5ad(DATA_DIR / DATASET / "processed" / f"adata_run_celldancer_{fold}.h5ad")
    del adata
    gc.collect()

def run_celldancer_time(fold,n_paths):
    adata = sc.read_h5ad(DATA_DIR / DATASET / "processed" / f"adata_run_celldancer_{fold}.h5ad")
    cellDancer_df = pd.read_csv(DATA_DIR / DATASET / "processed" / f"cellDancer_df_{fold}.csv")
    dt = 0.001
    t_total = {dt: 10000}
    n_repeats = 10
    cellDancer_df = cd.pseudo_time(cellDancer_df=cellDancer_df,
                                              grid=(30, 30),
                                              dt=dt,
                                              t_total=t_total[dt],
                                              n_repeats=n_repeats,
                                              speed_up=(60,60),
                                              n_paths = n_paths,
                                              psrng_seeds_diffusion=[i for i in range(n_repeats)],
                                              n_jobs=8)
    adata = export_velocity_to_dynamo(cellDancer_df,adata)
    subdf = cellDancer_df.iloc[0:adata.shape[0], :]
    subdf.index = adata.obs.index
    adata.obs['celldancer_time'] = subdf['pseudotime']
    if SAVE_DATA:
        adata.write_h5ad(DATA_DIR / DATASET / "processed" / f"adata_run_celldancer_{fold}.h5ad")
    del adata  # Remove reference to AnnData object
    gc.collect()

if K_FOLD == 0:
    fold = 'full'
    run_celldancer(fold)
else:
    for fold in range(K_FOLD):
        run_celldancer(fold)
    for item in os.listdir("."):
        if item.startswith("cellDancer") and os.path.isdir(item):
            shutil.rmtree(item)
    for item in os.listdir("."):
        if item.startswith("CellDancer") and os.path.isdir(item):
            shutil.rmtree(item)
            
    candidates = [3, 2, 1]
    for fold in range(1,K_FOLD):
        success = False
        for n_paths in candidates:
            try:
                run_celldancer_time(fold, n_paths)
                success = True
                break
            except Exception as e:
                print(f"fold {fold}, n_paths={n_paths} error: {e}")
        if not success:
            print(f"⚠️ fold {fold} all n_paths failed, skip")