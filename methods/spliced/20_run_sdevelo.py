import argparse
import sdevelo as sv
import scvelo as scv
import scanpy as sc
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

args = sv.Config()
args.cuda_device = "cuda:0"
args.vis_type_col = CLUSTER_KEY
args.vis_key = VIS_EMB

def run_sdevelo(fold, args):
    adata = sc.read_h5ad(DATA_DIR / DATASET / "processed" / f"adata_preprocessed_{fold}.h5ad")
    model = sv.SDENN(args, adata)
    adata = model.train(args.nEpochs)
    adata.layers['sdevelo_velocity'] = adata.layers.pop('sde_velocity')
    adata.obs['sdevelo_time'] = adata.obs['latent_time']
    if SAVE_DATA:
        adata.write_h5ad(DATA_DIR / DATASET / "processed" / f"adata_run_sdevelo_{fold}.h5ad")
    del adata  # Remove reference to AnnData object
    gc.collect()

if K_FOLD == 0:
    fold = 'full'
    run_sdevelo(fold)
else:
    for fold in range(K_FOLD):
        run_sdevelo(fold, args)
    