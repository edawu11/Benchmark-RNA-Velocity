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
  [26_run_multivelo.py]=py310_pt280
  [27_run_latentvelo_atac.py]=py310_pt260
  [28_run_graphvelo_std.py]=py310_pt212
)

# Ordered list of scripts to execute
ordered_scripts=(
  26_run_multivelo.py
  27_run_latentvelo_atac.py
  28_run_graphvelo_std.py
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
