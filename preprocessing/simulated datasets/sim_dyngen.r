#!/usr/bin/env Rscript

# ==============================================================================
# 1. Environment Configuration & Library Loading
# ==============================================================================
# Set CXX compiler path (Specific to your HPC environment)
Sys.setenv(CXX = "/share/home/kongch/anaconda3/envs/r433/bin/g++")

library(reticulate)

# Configure Python environment
# Ensure this path matches your specific conda environment
use_python("/share/home/kongch/anaconda3/envs/r433/bin/python", required = TRUE)

# Load required libraries
library(tidyverse)
library(dyngen)
library(anndata)
library(Rcpp)

# Set seed for reproducibility
seed <- 0
set.seed(seed)

# ==============================================================================
# 2. Define Simulation Parameters
# ==============================================================================
message("Configuring model parameters...")

# Create a bifurcating trajectory backbone
backbone <- backbone_bifurcating()

# Create a linear trajectory backbone
# backbone <- backbone_linear()

# Initialize the model configuration
config <- initialise_model(
  backbone = backbone,
  num_cells = 1500,
  num_tfs = nrow(backbone$module_info),
  num_targets = 100,
  num_hks = 100,
  verbose = FALSE,
  simulation_params = simulation_default(
    total_time = 1000,
    census_interval = 1, 
    ssa_algorithm = ssa_etl(tau = 0.01),
    experiment_params = simulation_type_wild_type(num_simulations = 100)
  )
)

# ==============================================================================
# 3. Execute Simulation
# ==============================================================================
message("Starting simulation...")

# Generate Transcription Factor (TF) network
model <- generate_tf_network(config)

# Sample target genes and housekeeping genes
model <- generate_feature_network(model)

# Generate kinetics (splicing/degradation rates etc.)
model <- generate_kinetics(model)

# Simulate gold standard (ground truth)
model <- generate_gold_standard(model)

# Simulate cells
model <- generate_cells(model)

# Simulate experiment effects (counts generation)
model <- generate_experiment(model)

message("Simulation completed successfully.")

# ==============================================================================
# 4. Save Output
# ==============================================================================
# Define output path
output_dir <- "simulation_data/dyngen_bifurcating"
output_file <- file.path(output_dir, paste0("adata_raw_",seed,".h5ad"))

# Create directory if it doesn't exist
if (!dir.exists(output_dir)) {
  message(paste("Creating output directory:", output_dir))
  dir.create(output_dir, recursive = TRUE)
}

# Convert dyngen model to AnnData object
message("Converting to AnnData object...")
ad <- as_anndata(model)

# Write to .h5ad file
message(paste("Saving output to:", output_file))
ad$write_h5ad(output_file)

message("Done.")