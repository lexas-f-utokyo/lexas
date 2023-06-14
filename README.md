# LEXAS: Lifescience EXperiment seArch and Suggestion

## Introduction

This repository contains the source code used to develop LEXAS (https://lexas.f.u-tokyo.ac.jp).

LEXAS curates the descriptions of biomedical experiments and suggests genes
that could be analyzed and a specific experimental method for the next experiment.



## Installation

1. Clone the source code.

```
git clone https://github.com/lexas-f-utokyo/lexas.git
```

2. Download required data from [google drive repository](https://drive.google.com/file/d/15hQMmr4cCejZj5HR2q03ieu1Q6pjVwjI/view?usp=sharing) in the root directory and unzip the tar.gz files.

```
tar -zxvf Repository.tar.gz
```

Your final project directory structure should look like this:
.
├── lexas
│   ├── 'Feature collection.ipynb'
│   ├── Main.ipynb
│   ├── __pycache__
│   ├── data
│   ├── README.md
│   ├── articles
│   ├── lexas
│   └── model
└── Repository
    ├── Pretrained_model
    ├── biobert
    └── feature

## Dependencies
- ahocorapy >= 1.6.1
- en_core_sci_sm >= 0.4.0
- numpy >= 1.19.5
- pandas >= 1.1.5
- scipy >= 1.5.4
- scispacy >= 0.4.0
- shap >= 0.39.0
- torch >= 1.9.1
- tqdm >= 4.61.1
- transformers >=4.3.3
- xgboost >=1.4.2


## Acknowledgment

LEXAS relies on many information sources:

- [PubMed Central](https://www.ncbi.nlm.nih.gov/pmc/)
- [HGNC](https://www.genenames.org/)
- [GO](http://geneontology.org/)
- [MGI](http://www.informatics.jax.org/)
- [HPO](https://hpo.jax.org/app/)
- [OMIM](https://www.omim.org/)
- [Orphanet](https://www.orpha.net/)
- [HPA](https://www.proteinatlas.org/)
- [BioGRID](https://thebiogrid.org/)
- [DepMap](https://depmap.org/)
- [ENCODE](https://www.encodeproject.org/)

We gratefully acknowledge these resources.
