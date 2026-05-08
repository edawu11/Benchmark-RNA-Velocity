import math
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt
from tivelo.main import tivelo
import argparse
import numpy as np
import gc
from pathlib import Path
import shutil
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

def run_tivelo_simple(fold):
    adata = sc.read_h5ad(DATA_DIR / DATASET / "processed" / f"adata_preprocessed_{fold}.h5ad")
    # sc.tl.leiden(adata, resolution=0.5) ### 在原本cluster失效的情况下使用
    # CLUSTER_KEY = 'leiden'
    adata_ = tivelo(adata, CLUSTER_KEY, VIS_EMB, data_name=DATASET, save_folder="results", njobs=8, tree_gene=None,
                show_fig=False, filter_genes=True, constrain=False, loss_fun="mse", only_s=False,
                alpha_1=1, alpha_2=0.1, batch_size=1024, n_epochs=100, velocity_key="velocity",
                adjust_DTI=False, show_DTI=False, cluster_edges=None,
                measure_performance=False)
    shutil.rmtree("results")
    adata_.layers["tivelo_simple_velocity"] = adata_.layers.pop("velocity")

    
    if SAVE_DATA:
        adata_.write_h5ad(DATA_DIR / DATASET / "processed" / f"adata_run_tivelo_simple_{fold}.h5ad")
    del adata_
    gc.collect()

if K_FOLD == 0:
    fold = 'full'
    run_tivelo_simple(fold)
else:
    for fold in range(K_FOLD):
        run_tivelo_simple(fold)