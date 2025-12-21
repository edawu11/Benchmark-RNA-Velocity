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
    adata = sc.read_h5ad(DATA_DIR / DATASET / "processed" / f"adata_preprocessed_label_{fold}.h5ad")
    preprocessor = dyn.pp.Preprocessor(cell_cycle_score_enable=True)
    preprocessor.preprocess_adata(adata, recipe='monocle', tkey='label_time', experiment_type='one-shot')
    dyn.tl.dynamics(adata)
    dyn.tl.reduceDimension(adata)
    dyn.tl.cell_velocities(adata, calc_rnd_vel=True, transition_genes=adata.var_names,)
    dyn.vf.VectorField(adata, basis=VIS_EMB.split('_')[1],map_topography=True, M=50,pot_curl_div=True)
    adata.layers["dynamo_m2_velocity"] = adata.layers.pop("velocity_T")
    adata.obs['dynamo_m2_time'] = adata.obs[f"{VIS_EMB.split('_')[1]}_ddhodge_potential"]
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