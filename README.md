# LEXAS: Lifescience EXperiment seArch and Suggestion

## Introduction

This repository contains the source code used to develop LEXAS (https://lexas.f.u-tokyo.ac.jp).

LEXAS curates the descriptions of biomedical experiments and suggests genes
that could be analyzed and a specific experimental method for the next experiment.



## Installation

1. Clone the source code.

```
git clone https://github.com/lexas-f-utokyo/lexas-tmp.git
```

2. Download required data from [google drive repository](https://drive.google.com/file/d/15hQMmr4cCejZj5HR2q03ieu1Q6pjVwjI/view?usp=sharing) in the root directory and unzip the tar.gz files.

```
tar -zxvf Repository.tar.gz
```


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

## LEXAS Search

In this section, we describe how to get information about gene-related experiments from PMC articles.

### 0. Downloading PMC articles

We have prepared seven sample articles.
If you want to get information about the experiments from the articles of your interest, please follow the steps below.

1. Download PMC articles of your interest from the [PMC FTP service](https://ftp.ncbi.nlm.nih.gov/pub/pmc/) and save them in the articles/.
2. Generate a list of PMCIDs in PMCID_list.txt with one ID per line.

### 1. Extracting the result sections

1. Run "1.result-extract.py" to extract only the result section from the articles, except the figure legends.
2. The output file will be saved in data/result_sections.txt.

### 2. Sentence segmentation

1. Run "2.sentence_segmentation.py" to perform sentence segmentation with scispacy.
2. The output file will be saved in data/result_sections_segmentation.txt

### 3. Extracting the sentences describing gene-related experiments

Extract sentences containing at least one gene and one experiment method using
a manually created experiment list (data/experiment_list.txt) and a human gene term list provided by HGNC.

The genes and experiment methods will be masked with special tokens "[GENE]" and "[EXPE]" for following relation-extraction task.

1. Run "3.Sentence-extraction.py"
2. The output file will be saved in data/masked_sentences.txt

### 4. Relation extraction between gene names and experiment methods

Using a bio-BERT model trained to extract the relation between genes and experiment methods,
you can obtain information about gene-related experiments.

1. Edit "4.Relation-extraction.ipynb" and select whether you use "[cpu]" or "[cuda]".
2. Run "4.Relation-extraction.ipynb".
3. The output file will be saved in data/masked_sentences_bert.txt

If you want to train a new relation extraction model, run "Relation-extraction-train.ipynb".

## LEXAS Suggestion

We have already prepared two pre-trained models, one uses features from databases and knowledgebases
and the other uses only databases.

If you want to use the pre-trained models, start from step 8.

### 5. Numbering the experiments

1. Run "5.py" to number the experiment for training a model.

### 6. Collecting gene features

To generate feature vectors, you have to create two dictionaries.
One dictionary has the gene name as the keys and the associated features as the values. 
The other dictionary has the name of the features as the keys and the number of the feature as the values.

Exceptionally, the Word2vec and DepMap features require only one dictionary.
They require only one dictionary with the gene names as the keys,
word vectors and cancer dependencies as the values, respectively.

1. Run "6.feature.py" to generate feature dictionaries.

### 7. Training an XGBoost model that suggests next experiments.

1. Edit "7. XGBoost_train.ipynb" and select whether you use all information sources or only databases.
2. Run "7. XGBoost_train.ipynb"

### 8. Suggesting the genes that could be examined in the next experiment

1. Run "8. XGBoost_suggestion.ipynb"

## License

This source code may be used for non-commercial purposes including  

- Research by academic institutions

- Non-commercial research, including research conducted within commercial organizations

- Personal use, including blog posts.

If you want to use the source code for a commercial purpose, please contact the following email addresses.

Miho Sakao : sakao [at_mark] todaitlo.jp



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
