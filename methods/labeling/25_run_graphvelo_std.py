import scanpy as sc
import scvelo as scv
import cellrank as cr
import argparse
from graphvelo.utils import mack_score, adj_to_knn
from graphvelo.graph_velocity import GraphVelo
from graphvelo.plot import gene_score_histogram
import matplotlib.pyplot as plt
import dynamo as dyn
import warnings
warnings.simplefilter("ignore")
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

def run_graphvelo_std(fold):
    file_path = DATA_DIR / DATASET / "processed" / f"adata_run_scvelo_dyn_{fold}.h5ad"
    
    if file_path.exists():
        adata = sc.read_h5ad(file_path)
    else:
        adata = sc.read_h5ad(DATA_DIR / DATASET / "processed" / f"adata_preprocessed_{fold}.h5ad")
        scv.tl.recover_dynamics(adata, n_jobs=8)
        scv.tl.velocity(adata, mode="dynamical", vkey='scvelo_dyn_velocity')
        scv.tl.velocity_graph(adata,vkey = 'scvelo_dyn_velocity')
        scv.tl.latent_time(adata,vkey='scvelo_dyn_velocity')
    indices, _ = adj_to_knn(adata.obsp['connectivities'])
    adata.uns['neighbors']['indices'] = indices
    mack_score(adata, ekey='Ms', vkey='scvelo_dyn_velocity', tkey='latent_time')
    confident_genes = adata.var['mack_score'].sort_values(ascending=False)[:100].index
    gv = GraphVelo(adata, gene_subset=confident_genes, xkey='Ms', vkey='scvelo_dyn_velocity')
    gv.train()
    adata.layers['graphvelo_std_velocity'] = gv.project_velocity(adata.layers['Ms'])
    adata.obsm[f"gv_{VIS_EMB.split('_')[1]}"] = gv.project_velocity(adata.obsm[VIS_EMB])
    adata.obsm[f"velocity_{VIS_EMB.split('_')[1]}"] = adata.obsm[f"gv_{VIS_EMB.split('_')[1]}"].copy()
    dyn.vf.VectorField(adata, basis=VIS_EMB.split('_')[1], M=100, pot_curl_div=True)
    adata.obs['graphvelo_std_time'] = adata.obs[f"{VIS_EMB.split('_')[1]}_ddhodge_potential"]
    if SAVE_DATA:
        adata.write_h5ad(DATA_DIR / DATASET / "processed" / f"adata_run_graphvelo_std_{fold}.h5ad")
    del adata  # Remove reference to AnnData object
    gc.collect()
if K_FOLD == 0:
    fold = 'full'
    run_graphvelo_std(fold)
else:
    for fold in range(K_FOLD):
        run_graphvelo_std(fold)