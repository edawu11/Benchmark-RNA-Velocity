import argparse
import anndata
import numpy as np
import scvelo as scv
import scanpy as sc
import sys
import torch
import os.path
import deepvelo as dv
import pickle as pickle
import matplotlib.pyplot as plt
import pandas as pd
from os.path import exists
import time
import gc

import torch # 如果pytorch安装成功即可导入
print(torch.cuda.is_available()) # 查看CUDA是否可用
print(torch.cuda.device_count()) # 查看可用的CUDA数量
print(torch.version.cuda) # 查看CUDA的版本

from pathlib import Path
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
    (MODEL_DIR / DATASET / METHOD).mkdir(parents=True, exist_ok=True)

def run_deepvelo(fold):
    adata = sc.read_h5ad(DATA_DIR / DATASET / "processed" / f"adata_preprocessed_{fold}.h5ad")
    trainer = dv.train(adata, dv.Constants.default_configs)
    adata.layers["deepvelo_velocity"] = adata.layers.pop("velocity")
    scv.tl.velocity_graph(adata, n_jobs=8)
    from deepvelo.utils.temporal import latent_time
    latent_time(adata)
    adata.obs['deepvelo_time'] = adata.obs['latent_time']
    if SAVE_DATA:
        adata.write_h5ad(DATA_DIR / DATASET / "processed" / f"adata_run_deepvelo_{fold}.h5ad")
    if torch.cuda.is_available():
        torch.cuda.empty_cache()  # Clear PyTorch's GPU memory cache
        torch.cuda.synchronize()
    del adata  # Remove reference to AnnData object
    gc.collect()

if K_FOLD == 0:
    fold = 'full'
    run_deepvelo(fold)
else:
    for fold in range(K_FOLD):
        run_deepvelo(fold)