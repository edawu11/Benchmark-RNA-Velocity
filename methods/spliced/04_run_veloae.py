import argparse
import os
import pickle
import numpy as np
import scvelo as scv
import scanpy as sc
import torch
import gc
from veloproj import *
SAVE_DATA = True

parser = get_parser()
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

args = parser.parse_args(args=['--lr',  '1e-6',
                               '--n-epochs', '20000', 
                               '--g-rep-dim', '100',
                               '--k-dim', '100',
                               '--data-dir', f'{DATA_DIR}/{DATASET}/processed/adata_run_scvelo_dyn_0.h5ad',
                               '--model-name', f'{DATASET}_0_model.cpt',
                               '--exp-name', f'CohAE_{DATASET}',
                               '--device', 'cuda:0',
                               '--gumbsoft_tau', '1',
                               '--nb_g_src', 'SU',
                               '--ld_nb_g_src', "X",
                               '--n_raw_gene', '2000',
                               '--n_conn_nb', '30',
                               '--n_nb_newadata', '30',
                               '--aux_weight', '1',
                               '--fit_offset_train', 'false',
                               '--fit_offset_pred', 'true',
                               '--use_offset_pred', 'false',
                               '--gnn_layer', 'GAT',
                               '--vis-key', VIS_EMB,
                               '--vis_type_col', CLUSTER_KEY,
                               '--scv_n_jobs', '10',
                               '--seed', str(SEED)
                              ])
torch.manual_seed(args.seed)
torch.cuda.manual_seed(args.seed)
np.random.seed(args.seed)
torch.backends.cudnn.deterministic = True
device = torch.device(args.device if args.device.startswith('cuda') and torch.cuda.is_available() else "cpu")

def main_AE(args, adata, tensor_vkey='scvelo_dyn_velocity'):
    spliced = adata.layers['Ms']  # moments of spliced
    unspliced = adata.layers['Mu'] # moments of unspliced
    tensor_s = torch.FloatTensor(spliced).to(device) 
    tensor_u = torch.FloatTensor(unspliced).to(device)
    tensor_x = torch.FloatTensor(adata.X.toarray()).to(device)
    tensor_v = torch.FloatTensor(adata.layers[tensor_vkey]).to(device)

    model = init_model(adata, args, device)

    inputs = [tensor_s, tensor_u]
    xyids = [0, 1]
    if args.use_x:
        inputs.append(tensor_x)

    model = fit_model(args, adata, model, inputs, tensor_v, xyids, device)
    return tensor_s, tensor_u, tensor_x  

def VeloAE_adata(adata):
    model = init_model(adata, args, device)
    model.load_state_dict(torch.load(args.model_name))
    model = model.to(device)
    model.eval()
    with torch.no_grad():
        
        spliced = adata.layers['Ms']
        unspliced = adata.layers['Mu']
        tensor_s = torch.FloatTensor(spliced).to(device)
        tensor_u = torch.FloatTensor(unspliced).to(device)
        tensor_x = torch.FloatTensor(adata.X.toarray()).to(device)
        

        x = model.encoder(tensor_x)
        s = model.encoder(tensor_s)
        u = model.encoder(tensor_u)
        
        v = estimate_ld_velocity(s, u, device=device, perc=[5, 95], 
                                 norm=args.use_norm, fit_offset=args.fit_offset_pred, 
                                 use_offset=args.use_offset_pred).cpu().numpy()
        x = x.cpu().numpy()
        s = s.cpu().numpy()
        u = u.cpu().numpy()
    
    adata = new_adata(adata, x, s, u, v, g_basis=args.ld_nb_g_src, 
                      n_nb_newadata=args.n_nb_newadata,new_v_key='veloae_velocity', X_emb_key=VIS_EMB)
    return adata

def run_veloae(fold, args):
    adata = sc.read_h5ad(args.data_dir)
    tensor_s, tensor_u, tensor_x = main_AE(args, adata)
    adata = VeloAE_adata(adata)
    if SAVE_DATA:
        adata.write_h5ad(DATA_DIR / DATASET / "processed" / f"adata_run_veloae_{fold}.h5ad")
    if torch.cuda.is_available():
        torch.cuda.empty_cache()  
        torch.cuda.synchronize()
    del adata
    gc.collect()
    
if K_FOLD == 0:
    fold = 'full'
    args.data_dir = DATA_DIR / DATASET / "processed" / f"adata_run_scvelo_dyn_{fold}.h5ad"
    args.model_name = DATA_DIR / DATASET / "processed" / f"{DATASET}_model_{fold}.cpt"
    args.exp_name = f"CohAE_{DATASET}_{fold}"
    run_veloae(fold, args)
else:
    for fold in range(K_FOLD):
        args.data_dir = DATA_DIR / DATASET / "processed" / f"adata_run_scvelo_dyn_{fold}.h5ad"
        args.model_name = DATA_DIR / DATASET / "processed" / f"{DATASET}_model_{fold}.cpt"
        args.exp_name = f"CohAE_{DATASET}_{fold}"
        run_veloae(fold, args)