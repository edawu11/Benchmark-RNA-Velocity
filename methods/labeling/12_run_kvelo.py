import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA #for creating PCAs
from sklearn.preprocessing import StandardScaler #for creating PCAs
import umap
from scipy.spatial import cKDTree
import sklearn as sk #used for L2 normalization
import time #to measure time of script running
import scanpy as sc
import gc
# import our own functions
import velocity
from scipy.sparse import issparse

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
METHOD = "kvelo"
SAVE_DATA = True
if SAVE_DATA:
    (DATA_DIR / DATASET / "processed").mkdir(parents=True, exist_ok=True)
    (MODEL_DIR / DATASET / METHOD).mkdir(parents=True, exist_ok=True)

def run_kvelo(fold):
    adata = sc.read_h5ad(DATA_DIR / DATASET / "processed" / f"adata_preprocessed_{fold}.h5ad")
    if issparse(adata.X):
        adata.X = adata.X.toarray()
    if issparse(adata.layers['spliced']):
        adata.layers['spliced'] = adata.layers['spliced'].toarray()
    if issparse(adata.layers['unspliced']):
        adata.layers['unspliced'] = adata.layers['unspliced'].toarray()
    if issparse(adata.layers['Ms']):
        adata.layers['Ms'] = adata.layers['Ms'].toarray()
    if issparse(adata.layers['Mu']):
        adata.layers['Mu'] = adata.layers['Mu'].toarray()
    adata.X = np.nan_to_num(adata.X, nan=0.0).astype(np.float32)
    adata.layers['spliced'] = np.nan_to_num(adata.layers['spliced'], nan=0.0).astype(np.float32)
    adata.layers['unspliced'] = np.nan_to_num(adata.layers['unspliced'], nan=0.0).astype(np.float32)
    adata.layers['Ms'] = np.nan_to_num(adata.layers['Ms'], nan=0.0).astype(np.float32)
    adata.layers['Mu'] = np.nan_to_num(adata.layers['Mu'], nan=0.0).astype(np.float32)
    minlim = 4
    us_genes = velocity.pp.filtering.get_high_us_genes(adata, minlim_u=minlim, minlim_s=minlim)
    us_genes = us_genes.astype(str)
    adata = adata[:,us_genes]
    velocity.pp.imputation.impute_counts(adata, n_neighbours = 30, layer_NN = 'spliced', n_pcs = 15)
    velocity.tl.fit.recover_reaction_rate_pars(adata, use_raw=False)
    likelihood_genes = adata.var['fit_likelihood'].sort_values(ascending=False)
    likelihood_genes = likelihood_genes.dropna()
    likelihood_genes = likelihood_genes[likelihood_genes>=0.1] ### the threshold need to be verify!!!!
    likelihood_genes = likelihood_genes.index.values
    likelihood_genes = likelihood_genes.astype(str)
    adata = adata[:, likelihood_genes]
    velocity.tl.fit.get_velocity(adata, use_raw=False, normalise=None)
    adata.layers["kvelo_velocity"] = adata.layers.pop("velocity")
    if SAVE_DATA:
        adata.write_h5ad(DATA_DIR / DATASET / "processed" / f"adata_run_kvelo_{fold}.h5ad")
    del adata  # Remove reference to AnnData object
    gc.collect()

if K_FOLD == 0:
    fold = 'full'
    run_kvelo(fold)
else:
    for fold in range(K_FOLD):
        run_kvelo(fold)
    