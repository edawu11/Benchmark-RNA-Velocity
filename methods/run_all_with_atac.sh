#!/bin/bash

# Exit on any error
set -e

# Load conda
source ~/miniconda3/etc/profile.d/conda.sh

# Create log folder
LOGDIR="./logs"
mkdir -p "$LOGDIR"

# Validate number of command-line arguments
if [ $# -lt 1 ]; then
    echo "Usage: $0 [dataset] [clusterkey] [X_embed] <K_FOLD> <seed>"
    exit 1
fi

# Assign command-line arguments to variables

dataset="$1"
clusterkey="${2:-}"
X_embed="${3:-}"
K_FOLD="${4:-}"
seed="${5:-}"

# Declare mapping: script => conda env
declare -A jobs=(
  [01_run_velocyto.py]=py310_pt212
  [02_run_scvelo_dyn.py]=py310_pt212
  [03_run_scvelo_stc.py]=py310_pt212
  [04_run_veloae.py]=py310_pt212
  [05_run_dynamo_m1.py]=py310_pt212
  [06_run_pyrovelocity_m1.py]=py311_pt221
  [07_run_pyrovelocity_m2.py]=py311_pt221
  [08_run_unitvelo_uni.py]=py310_tf213
  [09_run_unitvelo_ind.py]=py310_tf213
  [10_run_velovae_std.py]=py310_pt212
  [11_run_velovae_fullvb.py]=py310_pt212
  [12_run_kvelo.py]=py310_pt212
  [13_run_celldancer.py]=py37_pt110
  [14_run_velovi.py]=py310_pt260
  [15_run_latentvelo_std.py]=py310_pt260
  [16_run_sctour_mse.py]=py39_pt271
  [17_run_sctour_nb.py]=py39_pt271
  [18_run_sctour_zinb.py]=py39_pt271
  [19_run_deepvelo.py]=py39_pt113
  [20_run_sdevelo.py]=py310_pt212
  [21_run_svelvetvae.py]=py310_pt200
  [22_run_cell2fate.py]=py39_pt111
  [23_run_tivelo_std.py]=py39_pt250
  [24_run_tivelo_simple.py]=py39_pt250
  [25_run_graphvelo_std.py]=py310_pt212
  [26_run_multivelo.py]=py310_pt280
  [27_run_latentvelo_atac.py]=py310_pt260
  [28_run_graphvelo_std.py]=py310_pt212
)

# Ordered list of scripts to execute
ordered_scripts=(
  01_run_velocyto.py
  02_run_scvelo_dyn.py
  03_run_scvelo_stc.py
  04_run_veloae.py
  05_run_dynamo_m1.py
  06_run_pyrovelocity_m1.py
  07_run_pyrovelocity_m2.py
  08_run_unitvelo_uni.py
  09_run_unitvelo_ind.py
  10_run_velovae_std.py
  11_run_velovae_fullvb.py
  12_run_kvelo.py
  13_run_celldancer.py
  14_run_velovi.py
  15_run_latentvelo_std.py
  16_run_sctour_mse.py
  17_run_sctour_nb.py
  18_run_sctour_zinb.py
  19_run_deepvelo.py
  20_run_sdevelo.py
  21_run_svelvetvae.py
  22_run_cell2fate.py
  23_run_tivelo_std.py
  24_run_tivelo_simple.py
  25_run_graphvelo_std.py
)


# Log the configuration
echo "Configuration: dataset=$dataset, clusterkey=$clusterkey, X_embed=$X_embed, K_FOLD=$K_FOLD, seed=$seed"

# Loop through scripts in defined order
for script in "${ordered_scripts[@]}"; do
    env="${jobs[$script]}"
    log_file="$LOGDIR/${script%.py}.log"
    script_path="atac/$script"
    echo "Running $script in conda env: $env"
    source activate "$env"
    echo "=== Running $script ===" > "$log_file"
    echo "Config: dataset=$dataset, clusterkey=$clusterkey, X_embed=$X_embed, K_FOLD=$K_FOLD, seed=$seed" >> "$log_file"
    # Pass arguments, only including optional ones if non-empty
    cmd=(python "$script_path" --dataset "$dataset")
    [ -n "$clusterkey" ] && cmd+=(--clusterkey "$clusterkey")
    [ -n "$X_embed" ] && cmd+=(--X_embed "$X_embed")
    [ -n "$K_FOLD" ] && cmd+=(--K_FOLD "$K_FOLD")
    [ -n "$seed" ] && cmd+=(--seed "$seed")
    if ! "${cmd[@]}" >> "$log_file" 2>&1; then
        echo "[ERROR] $script failed. See $log_file for details."
        echo "[ERROR] Script failed: $script" >> "$log_file"
        # conda deactivate
        # exit 1
    fi
    conda deactivate
done

# Shutdown after all jobs finish
echo "All tasks complete. Shutting down system..."
# shutdown -h now   
