import numpy as np
from scipy import stats
from anndata import AnnData
from scipy.sparse import csr_matrix

def simulate_dropout(
        exp_counts: np.ndarray, 
        dropout_factor: float, 
        seed: int = 0) -> np.ndarray:
    """
    Simulate dropout using real expression count matrix.
    Args:
        exp_counts (np.ndarray): Expression count matrix with shape (n_cells, n_genes).
        dropout_factor (float): multiplier for global dropout rate control.
        seed (int): Random seed for reproducibility.
    """
    assert np.min(exp_counts) >= 0, "Expression counts must be non-negative."
    assert dropout_factor >=0 and dropout_factor <= 1, "Dropout factor must be in [0, 1]."
    np.random.seed(seed)
    num_cells = exp_counts.shape[0]
    
    # fill na values with 0
    exp_counts = exp_counts.astype(np.float32)
    exp_counts = np.nan_to_num(exp_counts, nan=0.0)
    
    # fit a logistic regression for parameter estimation
    gene_mean_expression = np.mean(exp_counts, axis=0)
    gene_wise_dropout_rate = np.sum(exp_counts == 0, axis=0) / num_cells
    log_gene_mean_expression = np.log(gene_mean_expression + 1e-12)

    eps = 1e-12
    gene_wise_dropout_rate = np.clip(gene_wise_dropout_rate, eps, 1 - eps)
    logit_dropout_rate = np.log(gene_wise_dropout_rate / (1 - gene_wise_dropout_rate))
    slope, intercept, *_ = stats.linregress(log_gene_mean_expression, logit_dropout_rate)

    k = slope
    x0 = -intercept / slope

    dropout_prob = 1 / (1 + np.exp(-k * (np.log(exp_counts + 1e-12) - x0)))
    dropout_prob = dropout_prob * dropout_factor

    dropout_mask = np.random.binomial(1, p=dropout_prob, size=exp_counts.shape)
    simulated_counts = exp_counts * (1 - dropout_mask)
    simulated_counts = np.nan_to_num(simulated_counts, nan=0.0).astype(np.float32)

    return simulated_counts, dropout_mask

def simulate_dropout_from_anndata(
        adata: AnnData,
        dropout_factor: float = 1,
        seed: int = 0
) -> AnnData:
    """ Simulate dropout on an AnnData object.

    Args:
        adata (AnnData): The AnnData object containing the expression data.
        dropout_factor (float): The factor by which to scale the dropout probability.
        seed (int): The random seed for reproducibility.

    Returns:
        AnnData: A new AnnData object with the simulated dropout applied.
    """
    assert "raw_spliced" in adata.layers, "The AnnData object must contain 'raw_spliced' layer."
    assert "raw_unspliced" in adata.layers, "The AnnData object must contain 'raw_unspliced' layer."

    raw_spliced = adata.layers["raw_spliced"].toarray().astype(np.float32)
    raw_unspliced = adata.layers["raw_unspliced"].toarray().astype(np.float32)
    spliced = adata.layers["spliced"].toarray().astype(np.float32)
    unspliced = adata.layers["unspliced"].toarray().astype(np.float32)

    simulated_raw_spliced, dropout_mask = simulate_dropout(
        exp_counts=raw_spliced,
        dropout_factor=dropout_factor,
        seed=seed
    )

    simulated_raw_unspliced = raw_unspliced * (1 - dropout_mask)
    simulated_raw_unspliced = np.nan_to_num(simulated_raw_unspliced, nan=0.0).astype(np.float32)

    simulated_spliced = spliced * (1 - dropout_mask)
    simulated_spliced = np.nan_to_num(simulated_spliced, nan=0.0).astype(np.float32)

    simulated_unspliced = unspliced * (1 - dropout_mask)
    simulated_unspliced = np.nan_to_num(simulated_unspliced, nan=0.0).astype(np.float32)

    adata_simulated = adata.copy()
    adata_simulated.layers["spliced"] = csr_matrix(simulated_spliced)
    adata_simulated.layers["unspliced"] = csr_matrix(simulated_unspliced)
    adata_simulated.layers["raw_spliced"] = csr_matrix(simulated_raw_spliced)
    adata_simulated.layers["raw_unspliced"] = csr_matrix(simulated_raw_unspliced)

    return adata_simulated


def downsample_gene_counts(
        exp_counts: np.ndarray,
        downsample_rate: float,
        seed: int = 0
) -> np.ndarray:
    """Downsample gene counts by a given rate.

    Args:
        exp_counts (np.ndarray): Expression count matrix with shape (n_cells, n_genes).
        downsample_rate (float): The rate at which to downsample the counts.
        seed (int): Random seed for reproducibility.
    Returns:
        np.ndarray: Downsampled expression count matrix.
    """
    assert downsample_rate >= 0 and downsample_rate <= 1, "Downsample rate must be in [0, 1]."
    np.random.seed(seed)
    exp_counts = np.nan_to_num(exp_counts, nan=0.0).astype(int)

    downsampled_counts = np.random.binomial(exp_counts, p=downsample_rate, size=exp_counts.shape)
    downsampled_counts = np.nan_to_num(downsampled_counts, nan=0.0).astype(np.float32)

    return downsampled_counts

def downsample_anndata(
        adata: AnnData,
        downsample_rate: float,
        seed: int = 0
) -> AnnData:
    """ Downsample gene counts in an AnnData object.    
    Args:
        adata (AnnData): The AnnData object containing the expression data.
        downsample_rate (float): The rate at which to downsample the counts.
        seed (int): The random seed for reproducibility.
    Returns:
        AnnData: A new AnnData object with the downsampled counts.
    """
    assert "raw_spliced" in adata.layers, "The AnnData object must contain 'raw_spliced' layer."
    assert "raw_unspliced" in adata.layers, "The AnnData object must contain 'raw_unspliced' layer."

    raw_spliced = adata.layers["raw_spliced"].toarray().astype(np.float32)
    raw_unspliced = adata.layers["raw_unspliced"].toarray().astype(np.float32)

    downsampled_raw_spliced = downsample_gene_counts(
        exp_counts=raw_spliced,
        downsample_rate=downsample_rate,
        seed=seed
    )

    downsampled_raw_unspliced = downsample_gene_counts(
        exp_counts=raw_unspliced,
        downsample_rate=downsample_rate,
        seed=seed
    )

    adata_downsampled = adata.copy()
    
    for key in adata.layers.keys():
        del adata_downsampled.layers[key]

    adata_downsampled.layers["spliced"] = csr_matrix(downsampled_raw_spliced)
    adata_downsampled.layers["unspliced"] = csr_matrix(downsampled_raw_unspliced)

    return adata_downsampled
