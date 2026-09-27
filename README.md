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
