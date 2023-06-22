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
- [Usage](#usage)
- [API](#api)
- [Contributing](#contributing)
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


## Usage

```
```

Note: The `license` badge image link at the top of this file should be updated with the correct `:user` and `:repo`.

### Any optional sections

## API

### Any optional sections

## More optional sections

## Contributing

See [the contributing file](CONTRIBUTING.md)!

PRs accepted.

Small note: If editing the Readme, please conform to the [standard-readme](https://github.com/RichardLitt/standard-readme) specification.

### Any optional sections

## License

[MIT © Richard McRichface.](../LICENSE)
