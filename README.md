# Comprehensive benchmarking of RNA velocity methods across single-cell datasets

## 📝 Description
We present a comprehensive benchmark of 19 computational RNA velocity tools covering 30 distinct methods. We systematically evaluate 25 splicing dynamics--based methods across eight evaluation tasks, designating directional consistency, temporal precision, negative control robustness, and sequencing depth stability as core tasks, while assessing five multimodal-enhanced methods specifically on the multimodal integration task. These assessments utilize 30 datasets encompassing 22 real-world and eight simulated scenarios.

## 🚀 RNA Velocity Methods
### ✂️ Splicing Dynamics–based 
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

### 🧬 Multimodal-enhanced 
| ID | Method | Modality | Paper | Code/Tutorial |
|:--:|:-------|:-------|:------|:--------------|
| 1 | protaccel | Protein abundance | [Protein velocity and acceleration from single-cell multiomics experiments](https://link.springer.com/article/10.1186/s13059-020-1945-3) | [PyPi](https://pypi.org/project/protaccel/0.301/); [Tutorial](https://github.com/pachterlab/GSP_2019) |
| 2 | MultiVelo | Chromatin accessibility | [Multi-omic single-cell velocity models epigenome–transcriptome interactions and improves cell fate prediction](https://www.nature.com/articles/s41587-022-01476-y) | [GitHub](https://github.com/welch-lab/MultiVelo); [Tutorial](https://multivelo.readthedocs.io/en/latest/) |
| 3 | PhyloVelo | Phylogenetic tree | [PhyloVelo enhances transcriptomic velocity field mapping using monotonically expressed genes](https://www.nature.com/articles/s41587-023-01887-5) | [GitHub](https://github.com/kunwang34/PhyloVelo); [Tutorial](https://phylovelo.readthedocs.io/) |
| 4 | VelvetVAE | Metabolic labeling | [Reconstructing developmental trajectories using latent dynamical systems and time-resolved transcriptomics](https://www.sciencedirect.com/science/article/pii/S2405471224001194?via%3Dihub) | [GitHub](https://github.com/rorymaizels/velvetVAE); [Tutorial](https://github.com/rorymaizels/Maizels2023aa/blob/main/analysis/A2.2_benchmarking/B01_velvet_benchmarking.py) |
| 5 | STT | Spatial transcriptomics | [Spatial transition tensor of single cells](https://www.nature.com/articles/s41592-024-02266-x) | [GitHub](https://github.com/cliffzhou92/STT/tree/release); [Tutorial](https://github.com/cliffzhou92/STT/tree/release/example_notebooks) |
| 6 | RegVelo | Gene regulatory network| [RegVelo: gene-regulatory-informed dynamics of single cells](https://www.biorxiv.org/content/10.1101/2024.12.11.627935v1) | [GitHub](https://github.com/theislab/RegVelo); [Tutorial](https://regvelo.readthedocs.io/en/latest/) |
| 7 | TFvelo | Transcription factors | [TFvelo: gene regulation inspired RNA velocity estimation](https://www.nature.com/articles/s41467-024-45661-w) | [GitHub](https://github.com/xiaoyeye/TFvelo) |
| 8 | spVelo | Spatial transcriptomics | [spVelo: RNA velocity inference for multi-batch spatial transcriptomics data](https://link.springer.com/article/10.1186/s13059-025-03701-8) | [GitHub](https://github.com/VivLon/spVelo); [Tutorial](https://github.com/VivLon/spVelo/blob/main/tutorial.ipynb) |
| 9 | TopoVelo | Spatial transcriptomics | [Topological velocity inference from spatial transcriptomic data](https://www.nature.com/articles/s41587-025-02688-8) | [GitHub](https://github.com/welch-lab/TopoVelo); [Tutorial](https://github.com/welch-lab/TopoVelo/tree/main/notebooks/tutorial) |

*Note: In our benchmark study, we specifically focus on methods enhanced by chromatin accessibility and metabolic labeling.*

## 📂 Dataset Information
| Data ID |  Dataset  | Access Method / Download Link |
|:-------:|:----------|:------------------------------|
| **Data 1** | Bone marrow | `scvelo.datasets.bonemarrow()` |
| **Data 2** | Cerebral cortex  | GEO: [GSE162170](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162170) |
| **Data 3** | Dentate gyrus | `scvelo.datasets.dentategyrus()` |
| **Data 4** | Gastrulation | `scvelo.datasets.gastrulation_erythroid()` |
| **Data 5** | Pancreas | `scvelo.datasets.pancreas()` |
| **Data 6** | scEU-seq organoid | `dynamo.sample_data.scEU_seq_organoid()` |
| **Data 7** | scNT-seq cortical neuron | `dynamo.scNT_seq_neuron_splicing()`<br>`dynamo.scNT_seq_neuron_labeling()` |
| **Data 8** | FUCCI U2OS | [Figshare Download](https://figshare.com/articles/dataset/FUCCI_U2OS_cells/22773761?file=40461716) |
| **Data 9** | Murine embryonic | GEO: [GSE142425](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142425) |
| **Data 10** | Primary visual cortex| GEO: [GSE190940](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190940) |
| **Data 11** | Reprogramming (Morris) | `cellrank.datasets.reprogramming_morris()` |
| **Data 12** | Reprogramming (Schiebinger)| GEO: [GSE122662](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122662) |
| **Data 13** | PBMC68k | `scvelo.datasets.pbmc68k()` |
| **Data 14** | PBMC194 | [10x Genomics Dataset](https://www.10xgenomics.com/datasets/pbm-cs-from-a-healthy-donor-whole-transcriptome-analysis-3-1-standard-4-0-0) |
| **Data 15** | PBMC381 | [10x Genomics Dataset](https://www.10xgenomics.com/datasets/8-k-pbm-cs-from-a-healthy-donor-2-standard-2-1-0) |
| **Data 16** | PBMC497 | [10x Genomics Dataset](https://www.10xgenomics.com/datasets/pbm-cs-from-a-healthy-donor-targeted-immunology-panel-3-1-standard-4-0-0) |
| **Data 17** | PBMC769 | [10x Genomics Dataset](https://www.10xgenomics.com/datasets/10-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0) |
| **Data 18** | ATAC brain | [10x Genomics Dataset](https://www.10xgenomics.com/datasets/fresh-embryonic-e-18-mouse-brain-5-k-1-standard-1-0-0)<br>Processed: [GitHub](https://github.com/welch-lab/MultiVelo/tree/main/Examples) |
| **Data 19** | ATAC HSPC | Processed: [Figshare (RNA)](https://doi.org/10.6084/m9.figshare.22575358.v1), [Figshare (ATAC)](https://doi.org/10.6084/m9.figshare.22575343.v1)<br>Raw: [GSE70677](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70677) |
| **Data 20** | ATAC skin | Processed: [Figshare (RNA)](https://doi.org/10.6084/m9.figshare.22575307.v1), [Figshare (ATAC)](https://doi.org/10.6084/m9.figshare.22575313.v1)<br>Raw: [GSE140203](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140203) |
| **Data 21** | FUCCI RPE1 | [Figshare Download](https://figshare.com/articles/dataset/FUCCI_RPE1_cells/22773776/1?file=40461722) |
| **Data 22** | scNT-seq hematopoiesis dynamics | `dyn.sample_data.hematopoiesis_raw()`|

*Note: Different datasets correspond to specific benchmark tasks; please refer to our benchmark paper for further details.*

## 🛠️ Implementation
For convenience, environment configuration files are provided in the `env` directory.

### ⚙️ Preprocessing
The `preprocessing` folder contains notebooks for:
* Preprocessing real-world datasets.
* Simulating sequencing depth variations.
* Processing quantification data.
* Generating synthetic data using dynamical models.

### ▶️ Running RNA Velocity Methods
We provide Python scripts for each RNA velocity method, paired with their corresponding Conda environments. These pipelines can be executed directly using the provided shell scripts.

For example, to run all splicing dynamics-based methods on **Data 1** (Bone Marrow), execute the following command:

```bash
bash run_all_spliced.sh bone_marrow clusters X_umap 3 1234
```

### ⚖️ Evaluation and Visualization
For evaluation and visualization, we developed the Python package [VeloEV](https://github.com/edawu11/VeloEV). Please refer to the [documentation](veloev.readthedocs.io/en/latest/) for detailed tutorials.

## 📊 Figure Reproduction
To facilitate reproducibility, the `Figure_reproduction` directory contains all summary results and notebooks required to regenerate the figures presented in the manuscript and supplementary materials.

## 📖 Reference

Yida Wu, Chuihan Kong, Xu Liao, Zhixiang Lin, Xiaobo Sun, Jin Liu. Comprehensive benchmarking of RNA velocity methods across single-cell datasets. Preprint. 2026.