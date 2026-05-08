import argparse
import os
os.environ["PYROVELOCITY_TESTING_FLAG"] = "False"
os.environ["PYROVELOCITY_LOG_LEVEL"] = "ERROR" 

import logging
logging.getLogger("jax").setLevel(logging.ERROR)
logging.getLogger("root").setLevel(logging.ERROR)

import yaml
from pyrovelocity.tasks.train import train_model
import scanpy as sc
import scvelo as scv
import gc
import numpy as np

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
MAX_EPOCHS = 1000
SUBSET_OBS = False
SUBSET_VARS = False

def compute_mean_vector_field(
    pos,
    adata,
    basis=VIS_EMB.split('_')[1],
    n_jobs=1,
    spliced="spliced_pyro",
    raw=False,
):
    sc.pp.neighbors(adata, use_rep="X_pca")

    adata.var["velocity_genes"] = True

    if spliced == "spliced_pyro":
        if raw:
            ut = pos["ut"]
            st = pos["st"]
            ut = ut / ut.sum(axis=-1, keepdims=True)
            st = st / st.sum(axis=-1, keepdims=True)
        else:
            ut = pos["ut"]
            st = pos["st"]
        adata.layers["spliced_pyro"] = st.mean(0).squeeze()
        # if ('u_scale' in pos) and ('s_scale' in pos): # TODO: two scale for Normal distribution
        if "u_scale" in pos:  # only one scale for Poisson distribution
            adata.layers["velocity_pyro"] = (
                ut * pos["beta"] / pos["u_scale"] - st * pos["gamma"]
            ).mean(0)
        else:
            if "beta_k" in pos:
                adata.layers["velocity_pyro"] = (
                    (ut * pos["beta_k"] - pos["st"] * pos["gamma_k"]).mean(0).squeeze()
                )
            else:
                adata.layers["velocity_pyro"] = (
                    ut * pos["beta"] - pos["st"] * pos["gamma"]
                ).mean(0)
        scv.tl.velocity_graph(
            adata, vkey="velocity_pyro", xkey="spliced_pyro", n_jobs=n_jobs
        )
    elif spliced in ["Ms"]:
        ut = adata.layers["Mu"]
        st = adata.layers["Ms"]
        if ("u_scale" in pos) and ("s_scale" in pos):
            adata.layers["velocity_pyro"] = (
                ut * pos["beta"] / (pos["u_scale"] / pos["s_scale"]) - st * pos["gamma"]
            ).mean(0)
        else:
            adata.layers["velocity_pyro"] = (
                ut * pos["beta"] - pos["st"] * pos["gamma"]
            ).mean(0)
        scv.tl.velocity_graph(adata, vkey="velocity_pyro", xkey="Ms", n_jobs=n_jobs)
    elif spliced in ["spliced"]:
        ut = adata.layers["unspliced"]
        st = adata.layers["spliced"]
        if ("u_scale" in pos) and ("s_scale" in pos):
            adata.layers["velocity_pyro"] = (
                ut * pos["beta"] / (pos["u_scale"] / pos["s_scale"]) - st * pos["gamma"]
            ).mean(0)
        else:
            adata.layers["velocity_pyro"] = (
                ut * pos["beta"] - pos["st"] * pos["gamma"]
            ).mean(0)
        scv.tl.velocity_graph(
            adata, vkey="velocity_pyro", xkey="spliced", n_jobs=n_jobs
        )

    scv.tl.velocity_embedding(adata, vkey="velocity_pyro", basis=basis)

def run_pyrovelocity_m1(fold):
    adata = sc.read_h5ad(DATA_DIR / DATASET / "processed" / f"adata_preprocessed_{fold}.h5ad")
    adata_model_pos = train_model(adata,
                              guide_type='auto_t0_constraint',
                              model_type='auto',
                              use_gpu="auto",
                              likelihood = "Poisson",
                              num_samples = 30,
                              log_every=100,
                              learning_rate = 0.01,
                              seed = SEED,
                              patient_improve=1e-3,
                              max_epochs=1000, 
                              offset=False,
                              library_size=True,
                              patient_init=45,
                              batch_size=4000, 
                              include_prior=True)
    compute_mean_vector_field(adata_model_pos[2], adata)
    adata.layers['pyrovelocity_m1_velocity'] = adata.layers.pop('velocity_pyro')
    cell_time = adata_model_pos[2]['cell_time']
    cell_time = np.squeeze(cell_time, axis = -1)
    cell_time_post_mean = np.mean(cell_time, axis = 0)
    adata.obs['pyrovelocity_m1_time'] = cell_time_post_mean
    if SAVE_DATA:
        adata.write_h5ad(DATA_DIR / DATASET / "processed" / f"adata_run_pyrovelocity_m1_{fold}.h5ad")
    del adata  # Remove reference to AnnData object
    gc.collect()

# run model 1
if K_FOLD == 0:
    fold = 'full'
    run_pyrovelocity_m1(fold)
else:
    for fold in range(K_FOLD):
        run_pyrovelocity_m1(fold)