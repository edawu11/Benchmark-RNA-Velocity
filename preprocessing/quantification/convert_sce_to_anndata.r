# ==============================================================================
# 1. Environment Configuration & Library Loading
# ==============================================================================
library(reticulate)
library(Matrix)
library(SingleCellExperiment)

# Configure Python environment
# Note: Update this path if running on a different machine
python_path <- "/root/autodl-tmp/conda/envs/scvelo_py39/python.exe"
use_python(python_path, required = TRUE)

# Import Python modules
anndata <- import("anndata")
np <- import("numpy")
scipy <- import("scipy.sparse")

# ==============================================================================
# 2. Helper Function: Convert R Sparse Matrix to Scipy CSR
# ==============================================================================
to_csr <- function(rmat) {
  stopifnot(inherits(rmat, "sparseMatrix"))
  
  # Convert Csparse -> Rsparse to access j/p slots directly
  if (inherits(rmat, "CsparseMatrix")) {
    rmat <- as(rmat, "dgRMatrix")
  }
  stopifnot(inherits(rmat, "RsparseMatrix"))
  
  # Extract data for CSR construction
  data    <- np$array(slot(rmat, "x"), dtype = "float64")       # non-zeros
  indices <- np$array(as.integer(slot(rmat, "j")), dtype="int32")  # column indices
  indptr  <- np$array(as.integer(slot(rmat, "p")), dtype="int32")  # row pointers
  dims    <- as.integer(slot(rmat, "Dim"))
  
  # Create Scipy CSR matrix
  scipy$csr_matrix(reticulate::tuple(data, indices, indptr),
                   shape = reticulate::tuple(dims[1], dims[2]))
}

# ==============================================================================
# 3. Data Loading
# ==============================================================================
# Define paths
base_path <- "D:/research/RNA_velocity/dataset/pancreas_quantification/"
file_name <- "Pancreas_sce_quantifications_combined.rds"
output_dir <- file.path(base_path, "raw")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

message("Loading SingleCellExperiment object...")
sce <- readRDS(file.path(base_path, file_name))

# Prepare common metadata (obs and var are shared across all outputs)
message("Preparing metadata...")
obs <- as.data.frame(colData(sce))
var <- data.frame(gene_id = rownames(sce), row.names = rownames(sce))

# ==============================================================================
# 4. Define Assays to Process
# ==============================================================================
# List of (spliced_name, unspliced_name, output_short_name)
assay_configs <- list(
  list(s = "dropest_spliced",  u = "dropest_unspliced",  name = "dropest"),
  list(s = "velocyto_spliced", u = "velocyto_unspliced", name = "velocyto"),
  list(s = "starsolo_spliced", u = "starsolo_unspliced", name = "starsolo"),
  list(s = "alevin_prepref_isoseparate_cdna_introns_gentrome_spliced",
       u = "alevin_prepref_isoseparate_cdna_introns_gentrome_unspliced", 
       name = "alevin"),
  list(s = "kallisto_bustools_prepref_isoseparate_exclude_spliced",
       u = "kallisto_bustools_prepref_isoseparate_exclude_unspliced", 
       name = "kallisto_bustools")
)

# ==============================================================================
# 5. Main Processing Loop
# ==============================================================================
for (config in assay_configs) {
  spliced_name <- config$s
  unspliced_name <- config$u
  method_name <- config$name
  
  message(paste0("Processing method: ", method_name, "..."))
  
  # Check if assays exist in the object
  if (!all(c(spliced_name, unspliced_name) %in% assayNames(sce))) {
    warning(paste("Assays not found for", method_name, "- skipping."))
    next
  }
  
  # 1. Extract and Transpose (Genes x Cells -> Cells x Genes for Python)
  # Note: scVelo/Scanpy expects Cells as rows (N x M)
  mat_spliced <- Matrix::t(assay(sce, spliced_name))
  mat_unspliced <- Matrix::t(assay(sce, unspliced_name))
  
  # 2. Ensure RsparseMatrix format
  if (!inherits(mat_spliced, "RsparseMatrix")) { mat_spliced <- as(mat_spliced, "RsparseMatrix") }
  if (!inherits(mat_unspliced, "RsparseMatrix")) { mat_unspliced <- as(mat_unspliced, "RsparseMatrix") }
  
  # 3. Convert to Python Scipy CSR
  py_spliced <- to_csr(mat_spliced)
  py_unspliced <- to_csr(mat_unspliced)
  
  # 4. Construct AnnData object
  adata <- anndata$AnnData(
    X = py_spliced,         # Use spliced counts as main X
    obs = obs,
    var = var,
    layers = dict(
      spliced = py_spliced,
      unspliced = py_unspliced
    )
  )
  
  # 5. Save to disk
  out_file <- file.path(output_dir, paste0("adata_", method_name, ".h5ad"))
  message(paste("Saving to:", out_file))
  adata$write_h5ad(out_file)
  
  # Clean up heavy objects to free memory
  rm(adata, py_spliced, py_unspliced, mat_spliced, mat_unspliced)
  gc()
}

message("All conversions completed.")