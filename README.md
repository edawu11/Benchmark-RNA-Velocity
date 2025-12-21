# Comprehensive benchmarking of RNA velocity methods across single-cell datasets

## Description
we systematically benchmark 19 RNA velocity tools (30 variants), including 25 splicing dynamic-based and 5 multimodal-enhanced methods, across 26 datasets. We ranked methods’ performance across four core tasks: directional consistency, temporal precision, negative control robustness, and sequencing depth stability. Additionally, we further characterized the methods by assessing quantification stability, multimodal integration, simulation experiment, and computational scalability. Overall, cellDancer, SDEvelo, and UniTVelo (uni) emerged as top methods in the core tasks, however, no single method performs optimally across all benchmark scenarios, highlighting task-specific strengths and weaknesses.

<!-- ![Pipeline](pipeline.png)

These methods utilize a variety of clustering strategies, including community detection-based, machine learning-based, deep learning-based, to analyze and cluster single-cell transcriptomic and proteomic data. Below is a collection of useful links and titles for each method. -->

## RNA velocity methods List
### Splicing dynamic-based 

|ID | Methods | Variants | Title |  |Code|
<!-- |:-------:|:------------------------------------------------------------------------:|---------------------------------------------------------------------------------------------------------------------|:--:| -->
|1| Velocyto | [RNA velocity of single cells](https://www.nature.com/articles/s41586-018-0414-6) |[Tutorial](https://velocyto.org/)|
|2| scVelo |[Generalizing RNA velocity to transient cell states through dynamical modeling](https://www.nature.com/articles/s41587-020-0591-3) | [GitHub](https://github.com/theislab/scvelo); [Tutorial](https://scvelo.readthedocs.io/en/stable/index.html)|
|3| VeloAE | [Representation learning of RNA velocity reveals robust cell transitions](https://www.pnas.org/doi/full/10.1073/pnas.2105859118) |[GitHub](https://github.com/qiaochen/VeloAE); [Tutorial](https://github.com/qiaochen/VeloAE/tree/main/notebooks)|
|4| Dynamo |  [Mapping transcriptomic vector fields of single cells](https://www.cell.com/cell/fulltext/S0092-8674(21)01577-4?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421015774%3Fshowall%3Dtrue) | [GitHub](https://github.com/aristoteleo/dynamo-release); [Tutorial](https://dynamo-release.readthedocs.io/en/latest/) |
|5| Pyro-Velocity  | [Pyro-Velocity: Probabilistic RNA Velocity inference from single-cell data](https://www.biorxiv.org/content/10.1101/2022.09.12.507691v2)  |[GitHub](https://github.com/pinellolab/pyrovelocity); [Tutorial](https://pyrovelocity.readthedocs.io/en/latest/)|
|6| UniTVelo      | [UniTVelo: temporally unified RNA velocity reinforces single-cell trajectory inference](https://www.nature.com/articles/s41467-022-34188-7) |[GitHub](https://github.com/StatBiomed/UniTVelo); [Tutorial](https://unitvelo.readthedocs.io/en/latest/)|



|7| [VeloVAE](https://github.com/campbio/celda)                            | Celda: a Bayesian model to perform co-clustering of genes into modules and cells into subpopulations using single-cell RNA-seq data. |[36110899](https://pubmed.ncbi.nlm.nih.gov/36110899/)|
|8| [$\kappa$-velo](https://bitbucket.org/scLCA/single_cell_lca/src/master/)     | Latent cellular analysis robustly reveals subtle diversity in large-scale single-cell RNA-seq data.                 |[31566233](https://pubmed.ncbi.nlm.nih.gov/31566233/)|
|9| [cellDancer](https://github.com/feiyoung/DR-SC.Analysis)                  | Joint dimension reduction and clustering analysis of single-cell RNA-seq and spatial transcriptomics data.          |[35349708](https://pubmed.ncbi.nlm.nih.gov/35349708/)|
|10| [veloVI](https://github.com/ttgump/scDCC)                             | Model-based deep embedding for constrained clustering analysis of single cell RNA-seq data.                        |[33767149](https://pubmed.ncbi.nlm.nih.gov/33767149/)|
|11| [LatentVelo](https://www.bioconductor.org/packages/release/bioc/html/SIMLR.html) | Visualization and analysis of single-cell RNA-seq data by kernel-based similarity learning.                 |[28263960](https://pubmed.ncbi.nlm.nih.gov/28263960/)|
|12| [scTour](https://scanpy.readthedocs.io/en/stable/)                    | From Louvain to Leiden: guaranteeing well-connected communities.                                                  |[30914743](https://pubmed.ncbi.nlm.nih.gov/30914743/)|
|13| [DeepVelo](https://github.com/zji90/TSCAN)                              | TSCAN: Pseudo-time reconstruction and evaluation in single-cell RNA-seq analysis.                                  |[27179027](https://pubmed.ncbi.nlm.nih.gov/27179027/)|
|14| [SDEvelo](https://github.com/shibiaowan/SHARP)                         | SHARP: hyperfast and accurate processing of single-cell RNA-seq data via ensemble random projection.               |[31992615](https://pubmed.ncbi.nlm.nih.gov/31992615/)|
|15| [cell2fate](https://github.com/juexinwang/scGNN)                         | scGNN is a novel graph neural network framework for single-cell RNA-Seq analyses.                                  |[33767197](https://pubmed.ncbi.nlm.nih.gov/33767197/)|
|16| [TIVelo](https://github.com/tinglabs/scAIDE)                         | scAIDE: clustering of large-scale single-cell RNA-seq data reveals putative and rare cell types.                   |[33575628](https://pubmed.ncbi.nlm.nih.gov/33575628/)|
|17| [GraphVelo](https://github.com/igrabski/sc-SHC)                         | Significance analysis for clustering with single-cell RNA-sequencing data.                                         |[37429993](https://pubmed.ncbi.nlm.nih.gov/37429993/)|

### Multimodal-enhanced 
|ID | Methods                                                                 | Title                                                                                                               |Code|
|:-------:|:------------------------------------------------------------------------:|---------------------------------------------------------------------------------------------------------------------|:--:|

|1| [protaccel](https://scanpy.readthedocs.io/en/stable/)                   | Fast unfolding of communities in large networks.  |[Link](https://iopscience.iop.org/article/10.1088/1742-5468/2008/10/P10008/meta)|
|2| [MultiVelo](https://github.com/jlakkis/CarDEC)                          | A joint deep learning model enables simultaneous batch effect correction, denoising, and clustering in single-cell transcriptomics. |[34035047](https://pubmed.ncbi.nlm.nih.gov/34035047/)|
|3| [PhyloVelo](https://github.com/biovault/SCHNELpy)                        | SCHNEL: scalable clustering of high dimensional single-cell data.                                                  |[33381821](https://pubmed.ncbi.nlm.nih.gov/33381821/)|
|4| [VelvetVAE](https://github.com/xuebaliang/scziDesk)                   | Deep soft K-means clustering with self-training for single-cell RNA sequence data.                                  |[33575592](https://pubmed.ncbi.nlm.nih.gov/33575592/)|
|5| [STT](https://www.bioconductor.org/packages/release/bioc/html/DepecheR.html) | Determination of essential phenotypic elements of clusters in high-dimensional entities—DEPECHE.        |[30845234](https://pubmed.ncbi.nlm.nih.gov/30845234/)|
|6| [RegVelo](https://github.com/SofieVG/FlowSOM)                        | FlowSOM: Using self-organizing maps for visualization and interpretation of cytometry data.                         |[25573116](https://pubmed.ncbi.nlm.nih.gov/25573116/)|
|7| [TFvelo](https://cran.r-project.org/web/packages/Spectrum/index.html) | Spectrum: fast density-aware spectral clustering for single and multi-omic data.                                 |[31501851](https://pubmed.ncbi.nlm.nih.gov/31501851/)|
|8| [spVelo](https://cole-trapnell-lab.github.io/monocle3/docs/clustering/) | The single-cell transcriptional landscape of mammalian organogenesis.                                          |[30787437](https://pubmed.ncbi.nlm.nih.gov/30787437/)|
|9| [TopoVelo](https://github.com/ZhenyiWangTHU/MarkovHC)  | MarkovHC: Markov hierarchical clustering for the topological structure of high-dimensional single-cell omics data with transition pathway and critical point detection. |[34850940](https://pubmed.ncbi.nlm.nih.gov/34850940/)|


## Implementation
To ensure the environment is properly configured for the specific requirements of each clustering method, begin by confirming that all dependencies and packages required by both Python and R scripts are installed. After the environment setup is complete, proceed to perform benchmarking tests for the clustering methods by running the following scripts:

### For Python-based Clustering Algorithms
- **Scripts**:
  - `./ClusteringAlgorithms_Python/script_rna_adt_clustering.py`
  - `./ClusteringAlgorithms_Python/script_integration_clustering.py`
- **Purpose**: These scripts implement and evaluate various Python-based clustering algorithms, providing insights into the efficiency and accuracy of each method.
- **Execution**: Run the scripts in your Python environment.
  ```bash
  python script_rna_adt_clustering.py
  python script_integration_clustering.py

### For R-based Clustering Algorithms
- **Scripts**:
  - `./ClusteringAlgorithms_R/script_rna_adt_clustering.R`
  - `./ClusteringAlgorithms_R/script_integration_clustering.R`
- **Purpose**: These scripts implement and evaluate various R-based clustering algorithms, providing insights into the efficiency and accuracy of each method.
- **Execution**: Run the scripts in your R environment.
  ```R
  source("script_rna_adt_clustering.R")
  source("script_integration_clustering.R")

## Acknowledgments
We would like to express our sincere appreciation to the developers of the single-cell clustering methods  and multi-omics integration approaches included in this benchmark. Their pioneering work has significantly advanced the field, making it possible to conduct in-depth analyses of multi-omics data.