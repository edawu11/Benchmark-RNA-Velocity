from pathlib import Path
import os
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import matplotlib.pyplot as plt
import argparse
import gc

DATA_DIR = Path("/root/autodl-tmp/dataset/ATAC")

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

mv.settings.VERBOSITY = 0

def run_multivelo(fold):
    adata = sc.read_h5ad(DATA_DIR / DATASET / "processed" / f"adata_preprocessed_{fold}.h5ad")
    adata_atac = sc.read_h5ad(DATA_DIR / DATASET / "processed" / f"adata_atac_preprocessed_{fold}.h5ad")
    adata_result = mv.recover_dynamics_chrom(adata,
                                         adata_atac,
                                         max_iter=5,
                                         init_mode="invert",
                                         parallel=True,
                                         n_jobs = 8,
                                         save_plot=False,
                                         rna_only=False,
                                         fit=True,
                                         n_anchors=500,
                                         extra_color_key=f"{CLUSTER_KEY}",
                                         embedding=f"{VIS_EMB}"
                                        )
    scv.tl.velocity_graph(adata_result)
    scv.tl.latent_time(adata_result)
    adata_result.layers['multivelo_velocity'] = adata_result.layers['velo_s'].copy()
    adata_result.obs['multivelo_time'] = adata_result.obs['latent_time'].copy()
    if SAVE_DATA:
        adata_result.write_h5ad(DATA_DIR / DATASET / "processed" / f"adata_run_multivelo_{fold}.h5ad")
    del adata_result, adata, adata_atac
    gc.collect()

if K_FOLD == 0:
    fold = 'full'
    run_multivelo(fold)
else:
    for fold in range(K_FOLD):
        run_multivelo(fold)
