# scDecipher: multi-view graph attention network for deciphering dominant cell communication assembly in tumor microenvironment
===========================================================================


[![license](https://img.shields.io/badge/python_-3.8.0_-blue)](https://www.python.org/)
[![license](https://img.shields.io/badge/torch_-1.12.0_-blue)](https://pytorch.org/)
[![license](https://img.shields.io/badge/scanpy_-1.9.0_-blue)](https://scanpy.readthedocs.io/en/stable/)
[![license](https://img.shields.io/badge/anndata_-0.8.0_-blue)](https://anndata-tutorials.readthedocs.io/en/latest/index.html/)
[![license](https://img.shields.io/badge/R_-4.2.2_-blue)](https://www.r-project.org/)

scDecipher is a toolkit to decipher the dominant cell communication assmebly (DCA) for downstream target genes based on scRNA-seq data by utilizing multi-view graph attention network. (i) scDecipher takes advantage of four state-of-the-art cell–cell communication analysis tools to systematically and reliably infer ligand–receptor (L-R) pairs from single-cell RNA sequencing data. (ii) Based on prior knowledge of L-R pairs and gene expression, scDecipher constructs a multi-view cell-cell communication network between different cell types with single-cell resolution by using an edge weighting strategy and filtering out edges with low specificity. (iii) scDecipher develops a multi-view graph attention network to predict the expression pattern of target genes or the functional status of receiver cells, and deciphers the dominant functional cell–cell communication by interpreting the trained model. The overview figure of scDecipher is shown as follows.


![Image text](https://github.com/jiboyalab/scDecipher/blob/main/IMG/overview.png)

## Table of Contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [Contributing](#contributing)
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
* Install [CellPhoneDB v3](https://github.com/ventolab/CellphoneDB) using ` pip install cellphonedb ` if you encounter any issue. 
* Install [CellChat v1.6.0](https://github.com/sqjin/CellChat/tree/master) using ` devtools::install_github("sqjin/CellChat") ` in the R environment if you encounter any issue.
* Install [NicheNet v1.1.0](https://github.com/saeyslab/nichenetr) using ` devtools::install_github("saeyslab/nichenetr") ` in the R environment if you encounter any issue.
* Install [ICELLNET](https://github.com/soumelis-lab/ICELLNET) using ` install_github("soumelis-lab/ICELLNET",ref="master", subdir="icellnet") ` in the R environment if you encounter any issue.


# Quick start
To reproduce our results:

1，infer ligand–receptor (L-R) pairs from single-cell RNA sequencing data
```
cellphonedb method statistical_analysis ./data/RCC_scRNA_P76_metadata.txt ./data/RCC_scRNA_P76_matrix.txt --counts-data=gene_name --threads 100 --output-path ./output/
```
**Arguments**:

| **Arguments** | **Detail** |
| --- | --- |
| **counts-data** | [ensembl or gene_name or hgnc_symbol] |
| **threads** | Max of threads to process the data. |
| **output-path** | Directory where the results will be allocated (the directory must exist). |

```
Rscript ./tools/run_cellchat.R --count ./data/RCC_scRNA_P76_matrix.txt --meta ./data/RCC_scRNA_P76_metadata.txt  --output ./output/

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

2，prioritize the dominant cell communication assmebly that regulates the target gene expression pattern
```

python ./src/process_final_lr.py --lr_cellphonedb ./output/process_cellphonedb_lr.csv --lr_cellchat ./output/process_cellchat_lr.csv --lr_nichenet ./output/process_nichenet_lr.csv --lr_icellnet ./output/process_icellchat_lr.csv --count ./data/RCC_scRNA_P76_matrix.txt --output ./output/final_lr.csv
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


## Contributing

Jiboya Xuliwen ..

## Contacts
If you have any questions or comments, please feel free to email: byj@hnu.edu.cn.

## License

[MIT © Richard McRichface.](../LICENSE)
