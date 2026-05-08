import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
from velovi import preprocess_data, VELOVI
import matplotlib.pyplot as plt
import seaborn as sns
import pickle as pickle
import matplotlib.pyplot as plt
import pandas as pd
import time
from os.path import exists
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

def add_velovi_outputs_to_adata(adata, vae):
    latent_time = vae.get_latent_time(n_samples=25)
    velocities = vae.get_velocity(n_samples=25, velo_statistic="mean")

    t = latent_time
    scaling = 20 / t.max(0)

    adata.layers["velovi_velocity"] = velocities / scaling
    adata.layers["latent_time_velovi"] = latent_time

    adata.var["fit_alpha"] = vae.get_rates()["alpha"] / scaling
    adata.var["fit_beta"] = vae.get_rates()["beta"] / scaling
    adata.var["fit_gamma"] = vae.get_rates()["gamma"] / scaling
    adata.var["fit_t_"] = (
        torch.nn.functional.softplus(vae.module.switch_time_unconstr)
        .detach()
        .cpu()
        .numpy()
    ) * scaling
    scaling = scaling.to_numpy() 
    adata.layers["fit_t"] = latent_time.values * scaling[np.newaxis, :]
    adata.var['fit_scaling'] = 1.0

def run_velovi(fold):
    adata = sc.read_h5ad(DATA_DIR / DATASET / "processed" / f"adata_preprocessed_{fold}.h5ad")
    VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
    vae = VELOVI(adata)
    vae.train()
    add_velovi_outputs_to_adata(adata, vae)
    scv.tl.latent_time(adata, min_likelihood=None)
    adata.obs['velovi_time'] = adata.obs['latent_time']
    if SAVE_DATA:
        adata.write_h5ad(DATA_DIR / DATASET / "processed" / f"adata_run_velovi_{fold}.h5ad")
    if torch.cuda.is_available():
        torch.cuda.empty_cache()  # Clear PyTorch's GPU memory cache
        torch.cuda.synchronize()
    del adata  # Remove reference to AnnData object
    gc.collect()

if K_FOLD == 0:
    fold = 'full'
    run_velovi(fold)
else:
    for fold in range(K_FOLD):
        run_velovi(fold)