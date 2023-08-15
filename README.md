# scDecipher: A multi-view graph attention network for deciphering dominant cell communication assembly with specific functional influence using single-cell RNA sequencing data
===========================================================================


[![license](https://img.shields.io/badge/python_-3.8.0_-blue)](https://www.python.org/)
[![license](https://img.shields.io/badge/torch_-1.12.0_-blue)](https://pytorch.org/)
[![license](https://img.shields.io/badge/scanpy_-1.9.0_-blue)](https://scanpy.readthedocs.io/en/stable/)
[![license](https://img.shields.io/badge/anndata_-0.8.0_-blue)](https://anndata-tutorials.readthedocs.io/en/latest/index.html/)
[![license](https://img.shields.io/badge/R_-4.2.2_-blue)](https://www.r-project.org/)

scDecipher is a toolkit to decipher the dominant cell communication assembly with specific functional influence by utilizing multi-view graph attention network. (i) scDecipher takes advantage of four state-of-the-art cell–cell communication analysis tools to systematically and reliably infer ligand–receptor (L-R) pairs from single-cell RNA sequencing data. (ii) Based on prior knowledge of L-R pairs and gene expression profile, scDecipher constructs a multi-view cell-cell communication network between different cell types at single-cell resolution by using an edge weighting strategy and filtering out edges with low specificity. (iii) scDecipher develops a multi-view graph attention network to predict the expression pattern of target genes or the functional status of receiver cells, and then deciphers the dominant cell communication assembly by interpreting the trained model. The overview figure of scDecipher is shown as follows.


![Image text](https://github.com/jiboyalab/scDecipher/blob/main/IMG/流程图_v2.svg)


The overview of scDecipher. **(a)** The schematic diagram of a full picture of cell-cell communication (CCC). LR: Ligand-Receptor; $E$: gene expression profile; $E_{0}$: baseline expression determined by cell type; $\Delta E$: expression change caused by CCC. $F$: functional states of malignant cells; $F_{0}$: baseline state determined by the cell itself; $\Delta F$: state change caused by CCC. **(b)** The construction of multi-view CCC network at single-cell resolution. First, four excellent CCC tools were applied to infer ligand-receptor pairs from single cell expression profiles. Second, an edge weighting strategy was applied to calculate the CCC strength at single-cell resolution based on the inferred LR pairs. Third, low specificity CCCs caused by widely expressed ligands (or receptors) or sequencing technology errors were filtered out. Fourth, the CCC graphs between each pair of different cell types were constructed, where the nodes in the graph represent single cells and the edges represent the communication strength between single cells calculated in the previous step. **(c)** The multi-view graph attention learning for predicting target gene expression or functional states of malignant cells. First, each CCC graph obtained from step b was trained by a graph convolutional network (GCN) and subjected to a nonlinear transformation by ReLU function, where the input contained the CCC strength matrix as the adjacency matrix and the one-hot encoding matrix obtained according to cell types as the initial feature matrix of the nodes. Second, the gene expressions or functional states obtained under different views were fused by an attention mechanism and then retrained using a multi-layer perceptron (MLP) to estimate $\Delta E/\Delta F$. Third, another MLP was applied to estimate  $E_{0}/F_{0}$ based on the cell-type matrix. Fourth, two estimations were summed as the final prediction results and trained iteratively to minimise the mean square error (MSE) with the true gene expression profiles or functional states. **(d)** The main function of scDecipher.
## Table of Contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [Contributing](#contributing)
- [Cite](#cite)
- [Contacts](#contacts)
- [License](#license)


## Installation

scDecipher is tested to work under:

```
* Python 3.8.0
* Torch 1.12.0
* Scanpy 1.9.0
* Anndata 0.8.0
* R 4.2.2
* Numpy 1.23.5
* Other basic python and r toolkits
```
### Installation of other dependencies 
**Notes:** These dependencies, in order to infer ligand–receptor (L-R) pairs from single-cell RNA sequencing data, can skip the installation process if you already have the LR result files (e.g. LR_P76.csv and LR_P915.csv provided in the data folder).
* Install [CellPhoneDB v3](https://github.com/ventolab/CellphoneDB) using ` pip install cellphonedb ` if you encounter any issue. 
* Install [CellChat v1.6.0](https://github.com/sqjin/CellChat/tree/master) using ` devtools::install_github("sqjin/CellChat") ` in the R environment if you encounter any issue.
* Install [NicheNet v1.1.0](https://github.com/saeyslab/nichenetr) using ` devtools::install_github("saeyslab/nichenetr") ` in the R environment if you encounter any issue.
* Install [ICELLNET](https://github.com/soumelis-lab/ICELLNET) using ` install_github("soumelis-lab/ICELLNET",ref="master", subdir="icellnet") ` in the R environment if you encounter any issue.


# Quick start
To reproduce our results:


**Notes:** We provide readers with 3 sets of data as detailed in the following data descriptions. Note that in order to reduce the computational overhead and to make it easier for readers to reproduce our code, we will use the smaller test data in the following tutorials. Processing of other single-cell RNA-Seq data follows the same core pipeline as the test data. Due to the large size of the data, we also uploaded them to the [Google Drive](https://drive.google.com/drive/folders/18sUpPlPuT9SBg-2rvurJ-IHjiKnhivI-?usp=drive_link).

## Data description

| File name  | Description |
| ------------- | ------------- |
| RCC_scRNA_P76_matrix.txt  | The single-cell gene expression matrix for patient 76 with advanced renal cell carcinoma. The origional data can be downloaded from the [Paper](https://singlecell.broadinstitute.org/single_cell/study/SCP1288/tumor-and-immune-reprogramming-during-immunotherapy-in-advanced-renal-cell-carcinoma?hiddenTraces=P55_scRNA%2CP90_scRNA%2CP906_scRNA%2CP912_scRNA%2CP913_scRNA%2CP916_scRNA%2CP76_scRNA#study-summary.) and the processed data by us can be obtained from the [Google Drive](https://drive.google.com/drive/folders/18sUpPlPuT9SBg-2rvurJ-IHjiKnhivI-?usp=drive_link)|
| RCC_scRNA_P76_metadata.txt  | The single-cell metadata, including cell type annotations, for patient 76 with advanced renal cell carcinoma. The origional data can be downloaded from the [Paper](https://singlecell.broadinstitute.org/single_cell/study/SCP1288/tumor-and-immune-reprogramming-during-immunotherapy-in-advanced-renal-cell-carcinoma?hiddenTraces=P55_scRNA%2CP90_scRNA%2CP906_scRNA%2CP912_scRNA%2CP913_scRNA%2CP916_scRNA%2CP76_scRNA#study-summary.) and the processed data by us can be obtained from the [Github/data](https://github.com/jiboyalab/scDecipher/tree/main/data) or [Google Drive](https://drive.google.com/drive/folders/18sUpPlPuT9SBg-2rvurJ-IHjiKnhivI-?usp=drive_link)|
| LR_P76.csv  | The integrated ligand-receptor results for patient 76 with advanced renal cell carcinoma obtained form 4 cell–cell communication analysis tools. The data can be obtained from the [Github/data](https://github.com/jiboyalab/scDecipher/tree/main/data) or [Google Drive](https://drive.google.com/drive/folders/18sUpPlPuT9SBg-2rvurJ-IHjiKnhivI-?usp=drive_link)|
| RCC_scRNA_P915_matrix.txt  | The single-cell gene expression matrix for patient 915 with advanced renal cell carcinoma. The origional data can be downloaded from the [Paper](https://singlecell.broadinstitute.org/single_cell/study/SCP1288/tumor-and-immune-reprogramming-during-immunotherapy-in-advanced-renal-cell-carcinoma?hiddenTraces=P55_scRNA%2CP90_scRNA%2CP906_scRNA%2CP912_scRNA%2CP913_scRNA%2CP916_scRNA%2CP76_scRNA#study-summary.) and the processed data by us can be obtained from the [Google Drive](https://drive.google.com/drive/folders/18sUpPlPuT9SBg-2rvurJ-IHjiKnhivI-?usp=drive_link)|
| RCC_scRNA_P915_metadata.txt  | The single-cell metadata, including cell type annotations, for patient 915 with advanced renal cell carcinoma. The origional data can be downloaded from the [Paper](https://singlecell.broadinstitute.org/single_cell/study/SCP1288/tumor-and-immune-reprogramming-during-immunotherapy-in-advanced-renal-cell-carcinoma?hiddenTraces=P55_scRNA%2CP90_scRNA%2CP906_scRNA%2CP912_scRNA%2CP913_scRNA%2CP916_scRNA%2CP76_scRNA#study-summary.) and the processed data by us can be obtained from the [Github/data](https://github.com/jiboyalab/scDecipher/tree/main/data) or [Google Drive](https://drive.google.com/drive/folders/18sUpPlPuT9SBg-2rvurJ-IHjiKnhivI-?usp=drive_link)|
| LR_P915.csv  | The integrated ligand-receptor results for patient 915 with advanced renal cell carcinoma obtained form 4 cell–cell communication analysis tools. The data can be obtained from the [Github/data](https://github.com/jiboyalab/scDecipher/tree/main/data) or [Google Drive](https://drive.google.com/drive/folders/18sUpPlPuT9SBg-2rvurJ-IHjiKnhivI-?usp=drive_link)|
| ScRNA_test_data_matrix.txt  | The single-cell gene expression matrix for test data to reduce the computational overhead and to make it easier for readers to reproduce our code. The origional data can be downloaded from the GSE175510 and the processed data by us can be obtained from the [Github/data](https://github.com/jiboyalab/scDecipher/tree/main/data) or [Google Drive](https://drive.google.com/drive/folders/18sUpPlPuT9SBg-2rvurJ-IHjiKnhivI-?usp=drive_link)|
| ScRNA_test_data_metadata.txt  | The single-cell metadata for test data including cell type annotations. The origional data can be downloaded from the GSE175510 and the processed data by us can be obtained from the [Github/data](https://github.com/jiboyalab/scDecipher/tree/main/data) or [Google Drive](https://drive.google.com/drive/folders/18sUpPlPuT9SBg-2rvurJ-IHjiKnhivI-?usp=drive_link)|
| LR_test_data.csv  | The integrated ligand-receptor results for test data. The data can be obtained from the [Github/data](https://github.com/jiboyalab/scDecipher/tree/main/data) or [Google Drive](https://drive.google.com/drive/folders/18sUpPlPuT9SBg-2rvurJ-IHjiKnhivI-?usp=drive_link)|
| P76_malignant_cell_states_gsva_mat.txt  | The activity scores calculated by gene set variation analysis (gsva) for 14 functional state in malignant cells of patient P76, which 14 functional state signatures of malignant cells were obtained from the [CancerSEA](http://biocc.hrbmu.edu.cn/CancerSEA/)|
| Test_data_malignant_cell_states_gsva_mat.txt  | The activity scores calculated by gene set variation analysis (gsva) for 14 functional state in malignant cells of test data, which 14 functional state signatures of malignant cells were obtained from the [CancerSEA](http://biocc.hrbmu.edu.cn/CancerSEA/)|

## 1，Infer ligand–receptor (L-R) pairs from single-cell RNA sequencing data
**Notes:** If you already have an LR result file or want to specify the LR yourself (e.g. LR_P76.csv, LR_P915.csv and LR_test_data.csv provided in the data folder), skip this step.
```
# The following program needs to be run in the cellphonedb environment, see [Cellphonedb](https://github.com/Teichlab/cellphonedb) for details on how to use it:
cellphonedb method statistical_analysis ./data/ScRNA_test_data_metadata.txt ./data/ScRNA_test_data_matrix.txt --counts-data=gene_name --iterations=10 --threads=100 --output-path=./output/
```
**Arguments**:

| **Arguments** | **Detail** |
| --- | --- |
| **counts-data** | [ensembl or gene_name or hgnc_symbol] Type of gene identifiers in the counts data |
| **iterations** | Number of iterations for the statistical analysis [1000] |
| **threads** | Number of threads to use. >=1 [4] |
| **output-path** | Directory where the results will be allocated (the directory must exist). |


```
[ ][CORE][15/08/23-10:17:35][INFO] Initializing SqlAlchemy CellPhoneDB Core
[ ][CORE][15/08/23-10:17:35][INFO] Using custom database at /home/jby2/.cpdb/releases/v4.0.0/cellphone.db
[ ][APP][15/08/23-10:17:35][INFO] Launching Method cpdb_statistical_analysis_local_method_launcher
[ ][APP][15/08/23-10:17:35][INFO] Launching Method _set_paths
[ ][APP][15/08/23-10:17:35][WARNING] Output directory (/home/jby2/HoloNet/github) exist and is not empty. Result can overwrite old results
[ ][APP][15/08/23-10:17:35][INFO] Launching Method _load_meta_counts
[ ][APP][15/08/23-10:17:37][INFO] Launching Method _check_counts_data
[ ][CORE][15/08/23-10:17:37][INFO] Launching Method cpdb_statistical_analysis_launcher
[ ][CORE][15/08/23-10:17:37][INFO] Launching Method _counts_validations
[ ][CORE][15/08/23-10:17:38][INFO] Launching Method get_interactions_genes_complex
[ ][CORE][15/08/23-10:17:38][INFO] [Cluster Statistical Analysis] Threshold:0.1 Iterations:10 Debug-seed:-1 Threads:100 Precision:3
[ ][CORE][15/08/23-10:17:39][INFO] Running Real Analysis
[ ][CORE][15/08/23-10:17:39][INFO] Running Statistical Analysis
[ ][CORE][15/08/23-10:17:43][INFO] Building Pvalues result
[ ][CORE][15/08/23-10:17:43][INFO] Building results
```


```
# The following program needs to be run in the R environment, see [CellChat](https://github.com/sqjin/CellChat), [NicheNet](https://github.com/saeyslab/nichenetr) and [ICELLNET](https://github.com/soumelis-lab/ICELLNET) for details on how to use it:

Rscript ./tools/run_cellchat.R --count ./data/ScRNA_test_data_matrix.txt --meta ./data/ScRNA_test_data_metadata.txt  --output ./output/
```
```
[1] "############ ------------- cellchat --------------- ############"
[1] ">>> loading library and data <<< [2023-08-15 10:48:39]"
[1] ">>> start CellChat workflow <<< [2023-08-15 10:48:50]"
[1] "Create a CellChat object from a data matrix"
The cell barcodes in 'meta' is  P2@CSF-0703-A1-1_GGTAATCA P2@CSF-0703-A1-1_CCTTCAAG P2@CSF-0703-A1-1_CACACTGA P2@CSF-0703-A1-2_GGTGGACT P2@CSF-0703-A1-2_TCAACGAC P2@CSF-0703-A2-1_GGTAATCA 
Set cell identities for the new CellChat object 
The cell groups used for CellChat analysis are  B CD8T Malignant Mono/Macro 
[1] ">>> Infer CCI network <<< [2023-08-15 10:49:05]"
triMean is used for calculating the average gene expression per cell group. 
[1] ">>> Run CellChat on sc/snRNA-seq data <<< [2023-08-15 10:49:05]"
  |======================================================================| 100%
[1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2023-08-15 10:49:46]"
[1] ">>> saving results <<< [2023-08-15 10:49:46]"
[1] ">>> done <<< [2023-08-15 10:49:48]"
Warning message:
In createCellChat(object = as.matrix(data.norm), group.by = "group",  :
  The cell barcodes in 'meta' is different from those in the used data matrix.
              We now simply assign the colnames in the data matrix to the rownames of 'mata'!
```


```
# The used ligand-target matrix, lr network and weighted networks of interacting cells can be downloaded from [Zenodo](https://zenodo.org/record/7074291).
Rscript ./tools/run_nichenet.R --count ./data/RCC_scRNA_P76_matrix.txt --meta ./data/RCC_scRNA_P76_metadata.txt  --output ./output/

Rscript ./tools/run_icellnet.R --count ./data/RCC_scRNA_P76_matrix.txt --meta ./data/RCC_scRNA_P76_metadata.txt  --output ./output/
```
**Arguments**:

| **Arguments** | **Detail** |
| --- | --- |
| **count** | Count matrix / normalized count matrix path. |
| **meta** | Meta data (celltypes annotation) path. |
| **output** | Directory where the results will be allocated. |

```
# Obtain the intersection of LR pairs output by 4 cellular communication tools, which are required to be found by at least 2 tools and have expression in scRNA-seq data.
python ./tools/process_final_lr.py --lr_cellphonedb ./output/process_cellphonedb_lr.csv --lr_cellchat ./output/process_cellchat_lr.csv --lr_nichenet ./output/process_nichenet_lr.csv --lr_icellnet ./output/process_icellchat_lr.csv --count ./data/RCC_scRNA_P76_matrix.txt --output ./output/final_lr.csv
```
**Arguments**:

| **Arguments** | **Detail** |
| --- | --- |
| **lr_cellphonedb** | The results of LR pairs output by cellphonedb. |
| **lr_cellchat** | The results of LR pairs output by cellchat. |
| **lr_nichenet** | The results of LR pairs output by nichenet. |
| **lr_icellnet** | The results of LR pairs output by icellnet. |
| **count** | Count matrix / normalized count matrix path. |
| **output** | The final results of LR pairs. |

## 2，prioritize the dominant cell communication assmebly that regulates the target gene expression pattern
```
python ./src/tutorials1/main.py --count ./data/RCC_scRNA_P76_matrix.txt --meta ./data/RCC_scRNA_P76_metadata.txt --lr_file ./output/final_lr.csv --gene CD8A --dca_rank_result ./output/CD8A_dca_rank_result.csv --ccc_ratio_result ./output/CD8A_ccc_ratio_result.csv
```
**Arguments**:

| **Arguments** | **Detail** |
| --- | --- |
| **count** | Count matrix / normalized count matrix path. |
| **meta** | Meta data (celltypes annotation) path. |
| **lr_file** | The final results of LR pairs. |
| **gene** | The specific target gene name  |
| **dca_rank_result** | The result of prioritize the dominant cell communication assmebly that regulates the target gene expression pattern. |
| **ccc_ratio_result** | The result of ratio of different cell types affected by cellular communication. |

Visualization of results:
<div align="center">
  <img src="https://github.com/jiboyalab/scDecipher/blob/main/IMG/cd8arank.png" alt="Editor" width="500">
</div>

===========================================================================

<div align="center">
  <img src="https://github.com/jiboyalab/scDecipher/blob/main/IMG/cd8adeltae.png" alt="Editor" width="400">
</div>

## 3，prioritize the dominant cell communication assmebly that regulates the key factors in specific cell type
```
python ./src/tutorials2/main.py --count ./data/RCC_scRNA_P76_matrix.txt --meta ./data/RCC_scRNA_P76_metadata.txt --lr_file ./output/final_lr.csv --gene FOLR2 --cell_type TAM --dca_rank_result ./output/FOLR2_dca_rank_result.csv --ccc_ratio_result ./output/FOLR2_ccc_ratio_result.csv
```
**Arguments**:

| **Arguments** | **Detail** |
| --- | --- |
| **count** | Count matrix / normalized count matrix path. |
| **meta** | Meta data (celltypes annotation) path. |
| **lr_file** | The final results of LR pairs. |
| **gene** | The specific target gene name.  |
| **cell_type** | The specific cell type (TAM:tumor-associated macrophages).  |
| **dca_rank_result** | The result of prioritize the dominant cell communication assmebly that regulates the target gene expression pattern. |
| **ccc_ratio_result** | The result of ratio of different cell types affected by cellular communication. |

Visualization of results:
<div align="center">
  <img src="https://github.com/jiboyalab/scDecipher/blob/main/IMG/folr2tam.png" alt="Editor" width="500">
</div>

## 4，prioritize the dominant cell communication assmebly that affected functional states of malignant cells
```
python ./src/tutorials3/main.py --count ./data/RCC_scRNA_P76_matrix.txt --meta ./data/RCC_scRNA_P76_metadata.txt --lr_file ./output/final_lr.csv --cell_type Malignant --cell_state EMT --dca_rank_result ./output/state_dca_rank_result.csv
```
**Arguments**:

| **Arguments** | **Detail** |
| --- | --- |
| **count** | Count matrix / normalized count matrix path. |
| **meta** | Meta data (celltypes annotation) path. |
| **lr_file** | The final results of LR pairs. |
| **cell_type** | The specific cell type.  |
| **cell_state** | [Angiogenesis; Apoptosis; CellCycle; Differentiation; DNAdamage; DNArepair; EMT; Hypoxia; Inflammation; Invasion; Metastasis; Proliferation; Quiescence; Stemness.]  |
| **dca_rank_result** | The result of prioritize the dominant cell communication assmebly that affected functional states of malignant cells. |


Visualization of results:
<div align="center">
  <img src="https://github.com/jiboyalab/scDecipher/blob/main/IMG/cellstate.png" alt="Editor" width="500">
</div>

## 5，clinical intervertion altered effect of cell communication on gene expression
```
python ./src/tutorials1/main.py --count ./data/RCC_scRNA_P915_matrix.txt --meta ./data/RCC_scRNA_P915_metadata.txt --lr_file ./output/final_lr.csv --gene CD8A --dca_rank_result ./output/P915_CD8A_dca_rank_result.csv --ccc_ratio_result ./output/P915_CD8A_ccc_ratio_result.csv
```
**Arguments**:

| **Arguments** | **Detail** |
| --- | --- |
| **count** | Count matrix / normalized count matrix path. |
| **meta** | Meta data (celltypes annotation) path. |
| **lr_file** | The final results of LR pairs. |
| **gene** | The specific target gene name  |
| **dca_rank_result** | The result of prioritize the dominant cell communication assmebly that regulates the target gene expression pattern. |
| **ccc_ratio_result** | The result of ratio of different cell types affected by cellular communication. |

Visualization of results:
<div align="center">
  <img src="https://github.com/jiboyalab/scDecipher/blob/main/IMG/cd8arankchange.png" alt="Editor" width="500">
</div>

===========================================================================





# Contributing

All authors were involved in the conceptualization of the proposed method. LWX and SLP conceived and supervised the project. BYJ and HTZ designed the study and developed the approach. BYJ and HTZ implemented and applied the method to microbial data. BYJ and HTZ analyzed the results. LWX and SLP contributed to the review of the manuscript before submission for publication. All authors read and approved the final manuscript.

# Cite
<p align="center">
  <a href="https://clustrmaps.com/site/1bpq2">
     <img width="200"  src="https://clustrmaps.com/map_v2.png?cl=ffffff&w=268&t=m&d=4hIDPHzBcvyZcFn8iDMpEM-PyYTzzqGtngzRP7_HkNs" />
   </a>
</p>

<p align="center">
  <a href="#">
     <img src="https://api.visitorbadge.io/api/visitors?path=https%3A%2F%2Fgithub.com%2Fjiboyalab%2FscDecipher&labelColor=%233499cc&countColor=%2370c168" />
   </a>
</p>


# Contacts
If you have any questions or comments, please feel free to email: byj@hnu.edu.cn.

# License

[MIT © Richard McRichface.](../LICENSE)
