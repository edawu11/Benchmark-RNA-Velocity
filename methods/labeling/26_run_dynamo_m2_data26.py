import argparse
import warnings
warnings.filterwarnings('ignore')
import dynamo as dyn
dyn.get_all_dependencies_version()
dyn.configuration.set_figure_params('dynamo', background='white')
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

def run_dynamo_m2(fold):
    adata = sc.read_h5ad(DATA_DIR / DATASET / "processed" / f"adata_preprocessed_{fold}.h5ad")
    del adata.layers['spliced'],adata.layers['unspliced'],adata.layers['raw_spliced'],adata.layers['raw_unspliced'],
    preprocessor = dyn.pp.Preprocessor(force_gene_list=adata.var_names)
    preprocessor.config_monocle_recipe(adata, n_top_genes=len(adata.var_names)) 
    preprocessor.preprocess_adata_monocle(adata,tkey='time', experiment_type='one-shot')
    dyn.tl.reduceDimension(adata)
    dyn.tl.moments(adata, group='time')
    adata.uns["pp"]["has_splicing"] = False
    dyn.tl.dynamics(adata,group='time',one_shot_method="sci_fate", model="deterministic")
    adata.layers["dynamo_m2_velocity"] = adata.layers.pop("velocity_T")
    dyn.tl.cell_velocities(
        adata,
        enforce=True,
        X=adata.layers["M_t"],
        V=adata.layers["dynamo_m2_velocity"],
        method="cosine",
    )
    if SAVE_DATA:
        adata.write_h5ad(DATA_DIR / DATASET / "processed" / f"adata_run_dynamo_m2_{fold}.h5ad")
    del adata  # Remove reference to AnnData object
    gc.collect()
if K_FOLD == 0:
    fold = 'full'
    run_dynamo_m2(fold)
else:
    for fold in range(K_FOLD):
        run_dynamo_m2(fold)