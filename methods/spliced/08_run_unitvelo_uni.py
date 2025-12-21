import argparse
import tensorflow as tf
print(tf.__version__)
print("GPU Available:", tf.config.list_physical_devices('GPU'))
from os.path import exists
import unitvelo as utv
import time
import pandas as pd
import numpy as np
import scanpy as sc
import scvelo as scv
import gc
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True 
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

# run unitvelo_unified
velo = utv.config.Configuration()
velo.R2_ADJUST = True
velo.IROOT = None
velo.FIT_OPTION = '1' # unified mode
velo.AGENES_R2 = 1

def run_unitvelo_uni(fold, velo):
    adata = sc.read_h5ad(DATA_DIR / DATASET / "processed" / f"adata_preprocessed_{fold}.h5ad")
    adata = utv.run_model(adata, CLUSTER_KEY, config_file=velo)
    adata.layers['unitvelo_uni_velocity'] = adata.layers.pop('velocity')
    adata.obs['unitvelo_uni_time'] = adata.obs['latent_time']
    if SAVE_DATA:
        adata.write_h5ad(DATA_DIR / DATASET / "processed" / f"adata_run_unitvelo_uni_{fold}.h5ad")
    del adata  # Remove reference to AnnData object
    gc.collect()

if K_FOLD == 0:
    fold = 'full'
    run_unitvelo_uni(fold, velo)
else:
    for fold in range(K_FOLD):
        run_unitvelo_uni(fold, velo)