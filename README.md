# TCGA-THCA miRNA Machine Learning Analysis

This repository contains the computational workflow used for the study:

**Machine Learning-Based Identification of a 30-miRNA Candidate Biomarker Panel for Thyroid Cancer Using TCGA-THCA Data**

## Overview

This study investigates miRNA expression profiles from the TCGA-THCA cohort using machine-learning approaches to identify candidate miRNA biomarkers capable of distinguishing primary thyroid tumor samples from solid tissue normal samples.

The analysis includes:

- TCGA-THCA miRNA expression data processing
- RPM normalization and log2(RPM+1) transformation
- Low-expression filtering
- Feature selection using ANOVA, mutual information (MI), and mRMR
- Comparison of original unbalanced data, SMOTE, and random undersampling
- Multiple machine-learning classifiers
- Patient-grouped nested cross-validation
- Performance evaluation using accuracy, F1-score, MCC, precision, recall, and ROC-AUC
- TOPSIS-based model ranking
- Selection of a 30-miRNA candidate biomarker panel

## Data source

miRNA expression data were obtained from the **TCGA-THCA** project through the National Cancer Institute Genomic Data Commons (GDC).

Raw TCGA data are not redistributed in this repository.

## Final model

The final candidate miRNA panel was derived from the highest-ranking hyperparameter-optimized configuration:

**Unbalanced data + mRMR feature selection + AdaBoost**

The model achieved:

- Accuracy: 0.991
- F1-score: 0.995
- MCC: 0.954
- Precision: 0.996
- Recall: 0.994
- ROC-AUC: 0.997
- TOPSIS score: 0.951
## Software environment

The computational environment corresponding to the analysis reported in the manuscript is summarized below.

### Python environment

- Python 3.9.23
- NumPy 2.0.2
- pandas 2.3.1
- scikit-learn 1.6.1
- SciPy 1.13.1
- imbalanced-learn 0.12.4
- matplotlib 3.9.4
- openpyxl 3.1.5
- mrmr-selection 0.2.8

The Python environment can be recreated using the `environment.yml` file provided in the root directory of this repository.

### R environment

- R 4.5.1
- TCGAbiolinks 2.36.0
- SummarizedExperiment 1.38.1
- dplyr 1.1.4
- stringr 1.5.2
- tibble 3.3.0
- readr 2.1.5
- purrr 1.1.0

The R script used for TCGA-THCA data acquisition and preprocessing is available in the `R/` directory.

## Repository structure

```text
R/          TCGA data acquisition and preprocessing
Python/     Machine-learning analysis
results/    Derived analysis results
Reproducibility
Software versions and the computational environment required to reproduce the final analysis will be provided in the accompanying environment file.
Data availability
The underlying TCGA-THCA data are publicly available from the NCI Genomic Data Commons.
License
License information will be added prior to publication.
