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

def run_dynamo_m1(fold):
    adata = sc.read_h5ad(DATA_DIR / DATASET / "processed" / f"adata_preprocessed_{fold}.h5ad")
    preprocessor = dyn.pp.Preprocessor(cell_cycle_score_enable=True)
    preprocessor.preprocess_adata(adata, recipe='monocle')
    dyn.tl.dynamics(adata, model='stochastic', est_method='negbin')
    dyn.tl.reduceDimension(adata, n_pca_components=30)
    dyn.tl.cell_velocities(adata, method='pearson', other_kernels_dict={'transform': 'sqrt'}, basis = VIS_EMB.split('_')[1])
    dyn.tl.cell_wise_confidence(adata)
    dyn.vf.VectorField(adata, basis=VIS_EMB.split('_')[1])
    dyn.ext.ddhodge(adata, basis=VIS_EMB.split('_')[1])
    adata.obs['dynamo_m1_time'] = adata.obs[f"{VIS_EMB.split('_')[1]}_ddhodge_potential"]
    adata.layers["dynamo_m1_velocity"] = adata.layers.pop("velocity_S")
    if SAVE_DATA:
        adata.write_h5ad(DATA_DIR / DATASET / "processed" / f"adata_run_dynamo_m1_{fold}.h5ad")
    del adata  # Remove reference to AnnData object
    gc.collect()
if K_FOLD == 0:
    fold = 'full'
    run_dynamo_m1(fold)
else:
    for fold in range(K_FOLD):
        run_dynamo_m1(fold)