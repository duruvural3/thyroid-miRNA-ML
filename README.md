# TCGA-THCA miRNA Machine Learning Analysis

This repository contains the computational workflows developed for the study:

**Machine Learning-Based Identification of a 30-miRNA Candidate Biomarker Panel for Thyroid Cancer Using TCGA-THCA Data**

## Overview

This study investigates miRNA expression profiles from the TCGA-THCA cohort using machine-learning approaches to identify candidate biomarkers capable of distinguishing primary thyroid tumor samples from solid tissue normal samples.

The study comprises two complementary computational approaches:

1. **miRNA Expression-Level Analysis:** Machine-learning analysis using TCGA-THCA miRNA expression quantification data.
2. **Mature miRNA Arm-Level (3p/5p) Analysis:** Machine-learning analysis of mature miRNAs with separate representation of 3p and 5p arms, using MIMAT-based isomiR aggregation.

Both approaches are integrated within the same research framework and are evaluated using feature-selection methods, machine-learning classifiers, patient-grouped nested cross-validation, and TOPSIS-based model ranking.

## Data Source

miRNA expression data were obtained from the **TCGA-THCA** project through the National Cancer Institute Genomic Data Commons (GDC).

The dataset comprises:

- 506 primary thyroid tumor samples
- 59 solid tissue normal samples
- 565 samples in total

Raw TCGA sequencing data are not redistributed in this repository.

## 1. miRNA Expression-Level Analysis

### Data Processing and Feature Selection

The computational workflow includes:

- TCGA-THCA miRNA expression data acquisition and preprocessing
- Reads-per-million (RPM) normalization
- log2(RPM+1) transformation
- Low-expression filtering
- Feature selection using analysis of variance (ANOVA), mutual information (MI), and minimum redundancy maximum relevance (mRMR)
- Comparison of unbalanced data, SMOTE oversampling, and random undersampling
- Evaluation of multiple machine-learning classifiers
- Patient-grouped nested cross-validation
- Performance evaluation using accuracy, F1-score, MCC, precision, recall, and ROC-AUC
- TOPSIS-based multi-criteria model ranking
- Identification of a 30-miRNA candidate biomarker panel

### Top-Ranked Model

The highest-ranking hyperparameter-optimized model configuration was:

**Unbalanced Data + mRMR Feature Selection + AdaBoost**

| Performance Metric | Value |
|---|---:|
| Accuracy | 0.991 |
| F1-score | 0.995 |
| MCC | 0.954 |
| Precision | 0.996 |
| Recall | 0.994 |
| ROC-AUC | 0.997 |
| TOPSIS score | 0.951 |

The corresponding computational scripts are available in the `R/` and `Python/` directories.

## 2. Mature miRNA Arm-Level (3p/5p) Analysis

This approach examines mature miRNA expression while distinguishing the 3p and 5p arms.

### Data Acquisition and Preprocessing

Mature miRNA expression data were retrieved from TCGA-THCA using GDC Isoform Expression Quantification files.

The R workflow performs the following operations:

1. Retrieves miRNA-Seq isoform expression data for primary tumor and solid tissue normal samples.
2. Retains isoforms annotated as mature miRNAs.
3. Extracts mature-miRNA MIMAT accession identifiers.
4. Aggregates read counts from isomiRs assigned to the same MIMAT identifier within each sample.
5. Maps MIMAT accessions to mature-miRNA names using miRBase version 21 annotations.
6. Aggregates counts at the mature-miRNA name level.
7. Retains one sample per patient and sample type, prioritizing the sample with the highest mature-miRNA library size.
8. Generates raw-count expression matrices for downstream computational analysis.

Tumor samples are labeled as 1 and normal samples as 0.

### Normalization and Feature Filtering

The Python workflow includes:

- RPM normalization based on mature-miRNA library sizes
- log2(RPM+1) transformation
- Low-expression filtering requiring RPM > 1 in at least 20% of samples

The mature-arm expression dataset comprises:

| Processing Stage | Number of Features |
|---|---:|
| Initial mature-miRNA features | 2,202 |
| Features retained after filtering | 533 |

### Machine-Learning Models

Ten machine-learning classifiers are evaluated:

- Support Vector Machine (SVM)
- Logistic Regression
- Extra Trees
- Linear Discriminant Analysis (LDA)
- Random Forest
- Gradient Boosting
- Stochastic Gradient Descent (SGD)
- K-Nearest Neighbors (KNN)
- Decision Tree
- AdaBoost

Three feature-selection methods are compared:

- ANOVA
- Mutual information (MI)
- Minimum redundancy maximum relevance (mRMR)

Three class-distribution strategies are evaluated:

- Unbalanced data
- SMOTE oversampling
- Random undersampling

Both default and hyperparameter-optimized model configurations are examined.

### Cross-Validation and Model Evaluation

Model evaluation is performed using patient-grouped nested cross-validation:

- Five outer folds
- Three inner folds for hyperparameter optimization

Classification performance is assessed using:

- Accuracy
- F1-score
- Matthews correlation coefficient (MCC)
- Precision
- Recall
- ROC-AUC

The Technique for Order Preference by Similarity to Ideal Solution (TOPSIS) is applied to rank model configurations using multiple performance criteria.

### Top-Ranked Model

The highest-ranking hyperparameter-optimized configuration in the mature-arm analysis was:

**SMOTE + mRMR Feature Selection + AdaBoost**

| Performance Metric | Value |
|---|---:|
| Accuracy | 0.986 |
| F1-score | 0.992 |
| MCC | 0.935 |
| Precision | 0.996 |
| Recall | 0.988 |
| ROC-AUC | 0.993 |
| TOPSIS score | 0.9933 |

The mature-arm workflow additionally produces feature-selection results, model comparison tables, biomarker direction information, and visualization outputs.

## 3. Software Environment

### Python Environment

The computational environment includes:

- Python 3.9.23
- NumPy 2.0.2
- pandas 2.3.1
- scikit-learn 1.6.1
- SciPy 1.13.1
- imbalanced-learn 0.12.4
- matplotlib 3.9.4
- openpyxl 3.1.5
- mrmr-selection 0.2.8

The Python environment is documented in the `environment.yml` file located in the repository root.

### R Environment

The R environment includes:

- R 4.5.1
- TCGAbiolinks 2.36.0
- SummarizedExperiment 1.38.1
- dplyr 1.1.4
- stringr 1.5.2
- tibble 3.3.0
- readr 2.1.5
- purrr 1.1.0

Additional R packages used in the mature-arm workflow include:

- tidyr
- AnnotationDbi
- miRBaseVersions.db

## 4. Repository Structure

```text
thyroid-miRNA-ML/
│
├── R/
│   └── miRNA expression data acquisition and preprocessing
│
├── Python/
│   └── miRNA expression-level machine-learning analysis
│
├── mature-arm-analysis/
│   ├── 01_TCGA_THCA_mature_arm_download.R
│   ├── 02_THCA_mature_arm_ML.ipynb
│   └── README.md
│
├── results/
│   └── Derived analysis results
│
├── environment.yml
├── CITATION.cff
├── LICENSE
└── README.md
```

## 5. Reproducibility

The computational workflows are provided to support methodological transparency and reproducibility.

### miRNA Expression-Level Workflow

The corresponding R and Python scripts are available in the `R/` and `Python/` directories.

### Mature miRNA Arm-Level Workflow

To reproduce the mature-arm analysis:

1. Execute `01_TCGA_THCA_mature_arm_download.R` to retrieve and preprocess TCGA-THCA mature-miRNA expression data.
2. Generate the mature-miRNA raw-count expression matrix.
3. Provide the generated expression matrix as input to `02_THCA_mature_arm_ML.ipynb`.
4. Execute the Python notebook to perform normalization, filtering, feature selection, model training, cross-validation, and TOPSIS ranking.
5. Examine the generated performance tables, selected features, and visualization outputs.

The computational environment is documented in the repository.

## 6. Data Availability

The underlying TCGA-THCA sequencing data are publicly available through the National Cancer Institute Genomic Data Commons:

https://portal.gdc.cancer.gov/

Raw sequencing data are not redistributed in this repository.

## 7. License

This repository is distributed under the terms specified in the `LICENSE` file.

## 8. Citation

Citation information and the associated repository DOI will be provided following the assignment of the corresponding research software release.
