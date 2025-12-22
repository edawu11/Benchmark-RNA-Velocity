# Comprehensive benchmarking of RNA velocity methods across single-cell datasets

## Description
<!-- we systematically benchmark 19 RNA velocity tools (30 variants), including 25 splicing dynamic-based and 5 multimodal-enhanced methods, across 26 datasets. We ranked methods’ performance across four core tasks: directional consistency, temporal precision, negative control robustness, and sequencing depth stability. Additionally, we further characterized the methods by assessing quantification stability, multimodal integration, simulation experiment, and computational scalability. Overall, cellDancer, SDEvelo, and UniTVelo (uni) emerged as top methods in the core tasks, however, no single method performs optimally across all benchmark scenarios, highlighting task-specific strengths and weaknesses. -->

<!-- ![Pipeline](pipeline.png)

These methods utilize a variety of clustering strategies, including community detection-based, machine learning-based, deep learning-based, to analyze and cluster single-cell transcriptomic and proteomic data. Below is a collection of useful links and titles for each method. -->

## RNA velocity methods List
### Splicing dynamics-based 
| ID | Method | Paper | Code/Tutorial |
|:--:|:-------|:------|:--------------|
| 1 | Velocyto | [RNA velocity of single cells](https://www.nature.com/articles/s41586-018-0414-6) | [Tutorial](https://velocyto.org/) |
| 2 | scVelo | [Generalizing RNA velocity to transient cell states through dynamical modeling](https://www.nature.com/articles/s41587-020-0591-3) | [GitHub](https://github.com/theislab/scvelo); [Tutorial](https://scvelo.readthedocs.io/en/stable/index.html) |
| 3 | VeloAE | [Representation learning of RNA velocity reveals robust cell transitions](https://www.pnas.org/doi/full/10.1073/pnas.2105859118) | [GitHub](https://github.com/qiaochen/VeloAE); [Tutorial](https://github.com/qiaochen/VeloAE/tree/main/notebooks) |
| 4 | Dynamo | [Mapping transcriptomic vector fields of single cells](https://www.cell.com/cell/fulltext/S0092-8674(21)01577-4?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421015774%3Fshowall%3Dtrue) | [GitHub](https://github.com/aristoteleo/dynamo-release); [Tutorial](https://dynamo-release.readthedocs.io/en/latest/) |
| 5 | Pyro-Velocity | [Pyro-Velocity: Probabilistic RNA Velocity inference from single-cell data](https://www.biorxiv.org/content/10.1101/2022.09.12.507691v2) | [GitHub](https://github.com/pinellolab/pyrovelocity); [Tutorial](https://pyrovelocity.readthedocs.io/en/latest/) |
| 6 | UniTVelo | [UniTVelo: temporally unified RNA velocity reinforces single-cell trajectory inference](https://www.nature.com/articles/s41467-022-34188-7) | [GitHub](https://github.com/StatBiomed/UniTVelo); [Tutorial](https://unitvelo.readthedocs.io/en/latest/) |
| 7 | VeloVAE | [Variational Mixtures of ODEs for Inferring Cellular Gene Expression Dynamics](https://arxiv.org/abs/2207.04166) | [GitHub](https://github.com/welch-lab/VeloVAE); [Tutorial](https://github.com/welch-lab/VeloVAE/blob/master/notebooks/velovae_example.ipynb) |
| 8 | $\kappa$-velo | [Towards reliable quantification of cell state velocities](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010031) | [GitHub](https://github.com/HaghverdiLab/velocity_package); [Tutorial](https://github.com/HaghverdiLab/velocity_notebooks/blob/main/tutorials/tutorial_02_kappavelo.ipynb) |
| 9 | cellDancer | [A relay velocity model infers cell-dependent RNA velocity](https://www.nature.com/articles/s41587-023-01728-5) | [GitHub](https://github.com/GuangyuWangLab2021/cellDancer); [Tutorial](https://guangyuwanglab2021.github.io/cellDancer_website/) |
| 10 | veloVI | [Deep generative modeling of transcriptional dynamics for RNA velocity analysis in single cells](https://www.nature.com/articles/s41592-023-01994-w) | [GitHub](https://github.com/YosefLab/velovi); [Tutorial](https://velovi.readthedocs.io/en/latest/tutorial.html) |
| 11 | LatentVelo | [Inferring single-cell transcriptomic dynamics with structured latent gene expression dynamics](https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(23)00225-4?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2667237523002254%3Fshowall%3Dtrue) | [GitHub](https://github.com/Spencerfar/LatentVelo); [Tutorial](https://github.com/Spencerfar/LatentVelo) |
| 12 | scTour | [scTour: a deep learning architecture for robust inference and accurate prediction of cellular dynamics](https://link.springer.com/article/10.1186/s13059-023-02988-9) | [GitHub](https://github.com/LiQian-XC/sctour/tree/main/sctour); [Tutorial](https://sctour.readthedocs.io/en/latest/) |
| 13 | DeepVelo | [DeepVelo: deep learning extends RNA velocity to multi-lineage systems with cell-specific kinetics](https://link.springer.com/article/10.1186/s13059-023-03148-9) | [GitHub](https://github.com/bowang-lab/DeepVelo); [Tutorial](https://github.com/bowang-lab/DeepVelo/blob/main/README.md) |
| 14 | SDEvelo | [Multivariate stochastic modeling for transcriptional dynamics with cell-specific latent time using SDEvelo](https://www.nature.com/articles/s41467-024-55146-5) | [GitHub](https://github.com/Liao-Xu/SDEvelo); [Tutorial](https://sdevelo.readthedocs.io/en/latest/tutorials/simulation.html) |
| 15 | cell2fate | [Cell2fate infers RNA velocity modules to improve cell fate prediction](https://www.nature.com/articles/s41592-025-02608-3) | [GitHub](https://github.com/BayraktarLab/cell2fate/tree/main); [Tutorial](https://cell2fate.readthedocs.io/en/latest/) |
| 16 | TIVelo | [TIVelo: RNA velocity estimation leveraging cluster-level trajectory inference](https://www.nature.com/articles/s41467-025-61628-x) | [GitHub](https://github.com/cuhklinlab/TIVelo); [Tutorial](https://tivelo.readthedocs.io/en/latest/) |
| 17 | GraphVelo | [GraphVelo allows for accurate inference of multimodal velocities and molecular mechanisms for single cells](https://www.nature.com/articles/s41467-025-62784-w) | [GitHub](https://github.com/xing-lab-pitt/GraphVelo); [Tutorial](https://graphvelo.readthedocs.io/en/latest/) |


### Multimodal-enhanced 
<!-- |ID | Methods                                                                 | Title                                                                                                               |Code|
|:-------:|:------------------------------------------------------------------------:|---------------------------------------------------------------------------------------------------------------------|:--:| -->

<!-- |1| [protaccel](https://scanpy.readthedocs.io/en/stable/)                   | Fast unfolding of communities in large networks.  |[Link](https://iopscience.iop.org/article/10.1088/1742-5468/2008/10/P10008/meta)|
|2| [MultiVelo](https://github.com/jlakkis/CarDEC)                          | A joint deep learning model enables simultaneous batch effect correction, denoising, and clustering in single-cell transcriptomics. |[34035047](https://pubmed.ncbi.nlm.nih.gov/34035047/)|
|3| [PhyloVelo](https://github.com/biovault/SCHNELpy)                        | SCHNEL: scalable clustering of high dimensional single-cell data.                                                  |[33381821](https://pubmed.ncbi.nlm.nih.gov/33381821/)|
|4| [VelvetVAE](https://github.com/xuebaliang/scziDesk)                   | Deep soft K-means clustering with self-training for single-cell RNA sequence data.                                  |[33575592](https://pubmed.ncbi.nlm.nih.gov/33575592/)|
|5| [STT](https://www.bioconductor.org/packages/release/bioc/html/DepecheR.html) | Determination of essential phenotypic elements of clusters in high-dimensional entities—DEPECHE.        |[30845234](https://pubmed.ncbi.nlm.nih.gov/30845234/)|
|6| [RegVelo](https://github.com/SofieVG/FlowSOM)                        | FlowSOM: Using self-organizing maps for visualization and interpretation of cytometry data.                         |[25573116](https://pubmed.ncbi.nlm.nih.gov/25573116/)|
|7| [TFvelo](https://cran.r-project.org/web/packages/Spectrum/index.html) | Spectrum: fast density-aware spectral clustering for single and multi-omic data.                                 |[31501851](https://pubmed.ncbi.nlm.nih.gov/31501851/)|
|8| [spVelo](https://cole-trapnell-lab.github.io/monocle3/docs/clustering/) | The single-cell transcriptional landscape of mammalian organogenesis.                                          |[30787437](https://pubmed.ncbi.nlm.nih.gov/30787437/)|
|9| [TopoVelo](https://github.com/ZhenyiWangTHU/MarkovHC)  | MarkovHC: Markov hierarchical clustering for the topological structure of high-dimensional single-cell omics data with transition pathway and critical point detection. |[34850940](https://pubmed.ncbi.nlm.nih.gov/34850940/)| -->


## Implementation
<!-- To ensure the environment is properly configured for the specific requirements of each clustering method, begin by confirming that all dependencies and packages required by both Python and R scripts are installed. After the environment setup is complete, proceed to perform benchmarking tests for the clustering methods by running the following scripts: -->

<!-- ### For Python-based Clustering Algorithms
- **Scripts**:
  - `./ClusteringAlgorithms_Python/script_rna_adt_clustering.py`
  - `./ClusteringAlgorithms_Python/script_integration_clustering.py`
- **Purpose**: These scripts implement and evaluate various Python-based clustering algorithms, providing insights into the efficiency and accuracy of each method.
- **Execution**: Run the scripts in your Python environment.
  ```bash
  python script_rna_adt_clustering.py
  python script_integration_clustering.py -->

<!-- ### For R-based Clustering Algorithms
- **Scripts**:
  - `./ClusteringAlgorithms_R/script_rna_adt_clustering.R`
  - `./ClusteringAlgorithms_R/script_integration_clustering.R`
- **Purpose**: These scripts implement and evaluate various R-based clustering algorithms, providing insights into the efficiency and accuracy of each method.
- **Execution**: Run the scripts in your R environment.
  ```R
  source("script_rna_adt_clustering.R")
  source("script_integration_clustering.R") -->

<!-- ## Acknowledgments
We would like to express our sincere appreciation to the developers of the single-cell clustering methods  and multi-omics integration approaches included in this benchmark. Their pioneering work has significantly advanced the field, making it possible to conduct in-depth analyses of multi-omics data. -->