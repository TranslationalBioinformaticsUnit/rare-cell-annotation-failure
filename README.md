# Rare-cell annotation failure: proximity and class imbalance

![R](https://img.shields.io/badge/R-%E2%89%A5%204.4.0-276DC3?logo=r&logoColor=white)
![License](https://img.shields.io/badge/license-MIT-green)

Analysis code for the manuscript **"Rare cell type annotation reveals geometry-associated vulnerabilities in single-cell foundation models"**.

A mixed-effects logistic regression predicts whether a rare cell type is poorly annotated (F1 < 0.5) from its **proximity** to the nearest common cell type in the learned representation and from **class imbalance**, with dataset as a random intercept. The script reports in-sample performance, feature contributions (SHAP) and leave-one-dataset-out cross-validation (LODO-CV).

## Contents

- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Reproducibility notes](#reproducibility-notes)
- [Data](#data)
- [Annotation frameworks benchmarked](#annotation-frameworks-benchmarked)
- [Datasets](#datasets)

## Requirements

R 4.4.0 or later and five packages. Versions used for the manuscript:

| Package | Version | Used for |
|---|---|---|
| [lme4](https://CRAN.R-project.org/package=lme4) | 1.1-36 | mixed-effects logistic regression (`glmer`) |
| [pROC](https://CRAN.R-project.org/package=pROC) | 1.18.5 | AUROC |
| [PRROC](https://CRAN.R-project.org/package=PRROC) | 1.4 | AUPRC |
| [fastshap](https://CRAN.R-project.org/package=fastshap) | 0.1.1 | Monte Carlo SHAP values |
| [dplyr](https://CRAN.R-project.org/package=dplyr) | 1.1.4 | data manipulation |

```r
install.packages(c("dplyr", "pROC", "PRROC", "lme4", "fastshap"))
# exact versions:
# remotes::install_version("lme4", "1.1-36")
```

## Quick start

```bash
git clone https://github.com/TranslationalBioinformaticsUnit/Rarecell-annotation-failure.git
cd Rarecell-annotation-failure
# Run analysis script
Rscript analysis.R
```

Run the script from the repository root. To analyze another representation, change `CASE` at the top of `analysis.R` (by default CellPLM_embeddings, options CellPLM_HVG, CELLama_embeddings, and CELLama_HVG).

## Reproducibility notes

- One seed (`SEED = 123`) is set at the top of the script and is used for SHAP sampling and the bootstrap, so keep the order of steps unchanged to reproduce the reported confidence interval.
- The final model is fitted with the default `glmer` optimizer; LODO-CV refits use `bobyqa` (`maxfun = 2e5`).
- LODO-CV pools out-of-fold predictions from all datasets because single datasets contain too few poorly annotated instances for a per-dataset AUROC.
- SHAP values are computed with `fastshap` on the predicted-probability scale, using fixed effects only.

## Data

`data/LR_data.csv` has one row per rare cell-type instance (a rare cell type within a dataset):

| Column | Description |
|---|---|
| `Case` | Analysis case, for example `CellPLM_embeddings` |
| `Dataset` | Dataset name (random intercept) |
| `low_f1` | `TRUE` if the cell type was poorly annotated (F1 < 0.5), otherwise `FALSE` |
| `Euclidean_min` | **Proximity**: Euclidean distance between the centroid of the rare cell type and that of the nearest common cell type; smaller means closer |
| `Imbalance_Ratio_Log_max` | **Class imbalance**: log-transformed imbalance ratio (common to rare cells), maximum over common cell types |

`data/Distances_rarecommon.csv` 
For every rare cell type, its distance to each common cell type of the same dataset, measured between cell-type centroids in a representation space (34,056 rows; four cases of 8,514 pairs each).

| Column | Description |
|---|---|
| `Pair` | `<rare cell type>vs<common cell type>` |
| `Dataset` | Dataset name |
| `Euclidean` | Euclidean distance between the two centroids |
| `Mahalanobis` | Mahalanobis (covariance-adjusted) distance between the two centroids |
| `Cosine_Similarity` | Cosine similarity between the two centroid vectors |
| `Imbalance_Ratio_Log` | Log-transformed imbalance ratio (cells of the common cell type divided by cells of the rare cell type) |
| `Case` | CellPLM_embeddings, CellPLM_HVG, CELLama_embeddings, and CELLama_HVG |
| `DatPair` | `Dataset_Pair` |
| `Rare_Celltype_ID` | `Dataset_<rare cell type>`, the key shared with `Benchmarking_results.csv` |


`data/Benchmarking_results.csv` 
F1-score of every cell type on the **held-out test partition** (30% of the cells of each dataset, stratified by cell type) for each annotation framework.

| Column | Description |
|---|---|
| `Tool` | Annotation framework (18 labels, see below) |
| `Dataset` | Dataset name |
| `Celltype` | Cell-type label; `macro avg` and `accuracy` rows are dataset-level summaries |
| `F1-score` | F1-score of the cell type on the test partition |
| `Id` | `Dataset_Celltype`, the key shared with `Distances_rarecommon.csv` |
| `Type` | `Rare` (< 1% of the cells of the dataset), `Common`, or `none` for summary rows |
| `SupportTest` | Number of cells of the cell type in the test partition |

The 18 `Tool` labels correspond to the 12 frameworks listed [below](#annotation-frameworks-benchmarked). UCE, scMulan and scPRINT, which produce embeddings, are evaluated with three downstream classifiers: logistic regression (`_LR`), 30-nearest neighbors (`_30KNN`) and a multilayer perceptron (`_MLP`). The benchmark summaries in the manuscript use the `_MLP`.


## Annotation frameworks benchmarked

Twelve models were benchmarked, each installed from its official GitHub repository at the commit listed below: nine foundation models, one transformer-based model and two machine-learning-based models.

| # | Model | Category | Repository | Commit |
|---|---|---|---|---|
| 1 | **CellLM** | Foundation model | [PharMolix/OpenBioMed](https://github.com/PharMolix/OpenBioMed) | [`557bcf9`](https://github.com/PharMolix/OpenBioMed/commit/557bcf9) |
| 2 | **CellPLM** | Foundation model | [OmicsML/CellPLM](https://github.com/OmicsML/CellPLM) | [`b59ee76`](https://github.com/OmicsML/CellPLM/commit/b59ee76) |
| 3 | **UCE (Universal Cell Embedding)** | Foundation model | [snap-stanford/UCE](https://github.com/snap-stanford/UCE) | [`7b31528`](https://github.com/snap-stanford/UCE/commit/7b31528) |
| 4 | **scMulan** | Foundation model | [SuperBianC/scMulan](https://github.com/SuperBianC/scMulan) | [`d64da5e`](https://github.com/SuperBianC/scMulan/commit/d64da5e) |
| 5 | **scPRINT** | Foundation model | [cantinilab/scPRINT](https://github.com/cantinilab/scPRINT) | [`2fcf5c7`](https://github.com/cantinilab/scPRINT/commit/2fcf5c7) |
| 6 | **SCimilarity** | Foundation model | [Genentech/scimilarity](https://github.com/Genentech/scimilarity) | [`4989ff8`](https://github.com/Genentech/scimilarity/commit/4989ff8) |
| 7 | **Nicheformer** | Foundation model | [theislab/nicheformer](https://github.com/theislab/nicheformer) | [`485cadb`](https://github.com/theislab/nicheformer/commit/485cadb) |
| 8 | **CellFM** | Foundation model | [biomed-AI/CellFM](https://github.com/biomed-AI/CellFM) | [`7655797`](https://github.com/biomed-AI/CellFM/commit/7655797) |
| 9 | **CELLama** | Foundation model | [portrai-io/CELLama](https://github.com/portrai-io/CELLama) | [`133c204`](https://github.com/portrai-io/CELLama/commit/133c204) |
| 10 | **TOSICA** | Transformer-based model | [JackieHanLab/TOSICA](https://github.com/JackieHanLab/TOSICA) | [`0aa2b32`](https://github.com/JackieHanLab/TOSICA/commit/0aa2b32) |
| 11 | **scMusketeers** | Machine-learning-based model | [AntoineCollin/scMusketeers](https://github.com/AntoineCollin/scMusketeers) | [`a806305`](https://github.com/AntoineCollin/scMusketeers/commit/a806305) |
| 12 | **scBalance** | Machine-learning-based model | [yuqcheng/scBalance](https://github.com/yuqcheng/scBalance) | [`5be134a`](https://github.com/yuqcheng/scBalance/commit/5be134a) |

## Datasets

The benchmark used 37 public scRNA-seq datasets.

| Dataset | Tissue | Source | Link |
|---|---|---|---|
| Covid | Blood, lung, saliva | CELLxGENE | [Collection](https://cellxgene.cziscience.com/collections/0a839c4b-10d0-4d64-9272-684c49a2c8ba) |
| Tabula | 72 tissues | CELLxGENE | [Collection](https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5) |
| Kumar | Breast | CELLxGENE | [Collection](https://cellxgene.cziscience.com/collections/4195ab4c-20bd-4cd3-8b3d-65601277e731) |
| Lymph | Lymph | CELLxGENE | [Collection](https://cellxgene.cziscience.com/collections/4429f25e-5c91-4981-9cf2-6a1044c36732) |
| Chen | Cornea | CELLxGENE | [Collection](https://cellxgene.cziscience.com/collections/0f7d022a-46c7-4e64-be4c-e34adbb78089) |
| HSPC | Blood | CELLxGENE | [Collection](https://cellxgene.cziscience.com/collections/5542eeb0-96ef-4ab9-95ea-eb6abc178461) |
| Heart atlas | Heart | CELLxGENE | [Collection](https://cellxgene.cziscience.com/collections/b52eb423-5d0d-4645-b217-e1c6d38b2e72) |
| Thymus atlas | Thymus | CELLxGENE | [Collection](https://cellxgene.cziscience.com/collections/fc19ae6c-d7c1-4dce-b703-62c5d52061b4) |
| Matuck oral | Oral and craniofacial atlas | CELLxGENE | [Collection](https://cellxgene.cziscience.com/collections/065ad318-59fd-4f8c-b4b1-66caa7665409) |
| Human retina | Retina | CELLxGENE | [Collection](https://cellxgene.cziscience.com/collections/5900dda8-2dc3-4770-b604-084eac1) |
| Spleen | Spleen | CELLxGENE | [Collection](https://cellxgene.cziscience.com/collections/854c0855-23ad-4362-8b77-6b1639e7a9fc) |
| Muscle | Intercostal muscle | CELLxGENE | [Collection](https://cellxgene.cziscience.com/collections/2d40e6a7-f2fd-49ba-9db9-6b97e4c6dad5) |
| Embryonic | Embryonic meninges | CELLxGENE | [Collection](https://cellxgene.cziscience.com/collections/7d66d871-091f-4602-9f42-85f86d2853e0) |
| Deng | T cell culture | CELLxGENE | [Collection](https://cellxgene.cziscience.com/collections/be21c2d1-2392-47d0-96fb-c625d115e0dc) |
| Zhang | Limb | CELLxGENE | [Collection](https://cellxgene.cziscience.com/collections/4fefa187-5d14-4f1e-915b-c892ed320aab) |
| Wong | Lung | CELLxGENE | [Collection](https://cellxgene.cziscience.com/collections/9f29fcd0-7075-4420-bd61-e2ec1aebad89) |
| Berg | Tracheal carina | CELLxGENE | [Collection](https://cellxgene.cziscience.com/collections/54004c5c-af08-4693-a606-73871b6ef989) |
| Barnes | Lung | CELLxGENE | [Collection](https://cellxgene.cziscience.com/collections/ec329aed-22bc-4d6e-8935-8282dcb1acac) |
| Zheng68K | Peripheral blood | 10x Genomics | [Dataset page](https://www.10xgenomics.com/datasets?menu%5Bproducts.name%5D=Single%20Cell%20Gene%20Expression) |
| Matuck dental | Dental | CELLxGENE | [Collection](https://cellxgene.cziscience.com/collections/065ad318-59fd-4f8c-b4b1-66caa7665409) |
| Bitzer | Kidney cortex | CELLxGENE | [Collection](https://cellxgene.cziscience.com/collections/a98b828a-622a-483a-80e0-15703678befd) |
| Human cortex | Dorsolateral prefrontal cortex | CELLxGENE | [Collection](https://cellxgene.cziscience.com/collections/ceb895f4-ff9f-403a-b7c3-187a9657ac2c) |
| Siletti | Cerebral cortex | CELLxGENE | [Collection](https://cellxgene.cziscience.com/collections/283d65eb-dd53-496d-adb7-7570c7caa443) |
| Matuck salivary | Minor salivary gland | CELLxGENE | [Collection](https://cellxgene.cziscience.com/collections/065ad318-59fd-4f8c-b4b1-66caa7665409) |
| Matuck buccal | Buccal mucosa | CELLxGENE | [Collection](https://cellxgene.cziscience.com/collections/065ad318-59fd-4f8c-b4b1-66caa7665409) |
| NeurIPS | Peripheral blood | Kaggle | [Competition data](https://www.kaggle.com/competitions/open-problems-multimodal/data) |
| Xpand | Bone marrow | EGA | [EGAS50000000750](https://ega-archive.org/studies/EGAS50000000750) |
| Deprez | Respiratory epithelium | Zenodo | [Record 12581306](https://zenodo.org/records/12581306) |
| MS | Cortical brain | EMBL-EBI | [E-HCAD-35](https://www.ebi.ac.uk/gxa/sc/experiments/E-HCAD-35) |
| Matuck gland | Gland atlas | CELLxGENE | [Collection](https://cellxgene.cziscience.com/collections/065ad318-59fd-4f8c-b4b1-66caa7665409) |
| hPancreas | Pancreas | GEO | [GSM2230757](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2230757) |
| PBMC | Peripheral blood | 10x Genomics | [Dataset page](https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets/1.0.0/pbmc_granulocyte_sorted_10k) |
| Matuck palate | Palate | CELLxGENE | [Collection](https://cellxgene.cziscience.com/collections/065ad318-59fd-4f8c-b4b1-66caa7665409) |
| Baron human | Pancreas | GEO | [GSE84133](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84133) |
| MacParland | Liver | GEO | [GSE115469](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115469) |
| Darmanis | Cerebral cortex | CELLxGENE | [Collection](https://cellxgene.cziscience.com/collections/558385a4-b7b7-4eca-af0c-9e54d010e8dc) |
| Muraro | Pancreas | CELLxGENE | [Collection](https://cellxgene.cziscience.com/collections/6e8c5415-302c-492a-a5f9-f29c57ff18fb) |

