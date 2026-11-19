TCGA-THCA miRNA Machine Learning Analysis

DOI

This repository contains the computational workflows developed for the study:

Comparative Identification of Locus-Level and Mature-Arm miRNA Candidate Signatures for Thyroid Carcinoma Using Patient-Grouped Machine Learning in TCGA-THCA

Overview

This study investigates miRNA expression profiles from the TCGA-THCA cohort using machine-learning approaches to identify candidate biomarkers capable of distinguishing primary thyroid tumor samples from solid tissue normal samples.

The study comprises two complementary computational approaches:

miRNA Expression-Level Analysis: Machine-learning analysis using TCGA-THCA miRNA expression quantification data.
Mature miRNA Arm-Level (3p/5p) Analysis: Machine-learning analysis of mature miRNAs with separate representation of 3p and 5p arms, using MIMAT-based isomiR aggregation.

Both approaches are integrated within the same research framework and are evaluated using feature-selection methods, machine-learning classifiers, patient-grouped nested cross-validation, and TOPSIS-based model ranking.

Model performance results, the derived candidate panels and their interpretation are reported in the associated manuscript; they are not duplicated here. Running the notebooks regenerates the corresponding performance tables, selected-feature sets and figures.

Data Source

miRNA expression data were obtained from the TCGA-THCA project through the National Cancer Institute Genomic Data Commons (GDC).

The dataset comprises:

506 primary thyroid tumor samples
59 solid tissue normal samples
565 samples in total

Raw TCGA sequencing data are not redistributed in this repository.

1. miRNA Expression-Level Analysis
Data Processing and Feature Selection

The computational workflow includes:

TCGA-THCA miRNA expression data acquisition and preprocessing
Reads-per-million (RPM) normalization
log2(RPM+1) transformation
Low-expression filtering
Feature selection using analysis of variance (ANOVA), mutual information (MI), and minimum redundancy maximum relevance (mRMR)
Comparison of unbalanced data, SMOTE oversampling, and random undersampling
Evaluation of multiple machine-learning classifiers
Patient-grouped nested cross-validation
Performance evaluation using accuracy, F1-score, MCC, precision, recall, and ROC-AUC
TOPSIS-based multi-criteria model ranking
Identification of a 30-miRNA candidate biomarker panel

The corresponding computational scripts are available in the R/ and Python/ directories.

2. Mature miRNA Arm-Level (3p/5p) Analysis

This approach examines mature miRNA expression while distinguishing the 3p and 5p arms.

Data Acquisition and Preprocessing

Mature miRNA expression data were retrieved from TCGA-THCA using GDC Isoform Expression Quantification files.

The R workflow performs the following operations:

Retrieves miRNA-Seq isoform expression data for primary tumor and solid tissue normal samples.
Retains isoforms annotated as mature miRNAs.
Extracts mature-miRNA MIMAT accession identifiers.
Aggregates read counts from isomiRs assigned to the same MIMAT identifier within each sample.
Maps MIMAT accessions to mature-miRNA names using miRBase version 21 annotations.
Aggregates counts at the mature-miRNA name level.
Retains one sample per patient and sample type, prioritizing the sample with the highest mature-miRNA library size.
Generates raw-count expression matrices for downstream computational analysis.

Tumor samples are labeled as 1 and normal samples as 0.

Normalization and Feature Filtering

The Python workflow includes:

RPM normalization based on mature-miRNA library sizes
log2(RPM+1) transformation
Low-expression filtering requiring RPM > 1 in at least 20% of samples

The mature-arm expression dataset comprises:

Processing Stage	Number of Features
Initial mature-miRNA features	2,202
Features retained after filtering	533
Machine-Learning Models

Ten machine-learning classifiers are evaluated:

Support Vector Machine (SVM)
Logistic Regression
Extra Trees
Linear Discriminant Analysis (LDA)
Random Forest
Gradient Boosting
Stochastic Gradient Descent (SGD)
K-Nearest Neighbors (KNN)
Decision Tree
AdaBoost

Three feature-selection methods are compared:

ANOVA
Mutual information (MI)
Minimum redundancy maximum relevance (mRMR)

Three class-distribution strategies are evaluated:

Unbalanced data
SMOTE oversampling
Random undersampling

Both default and hyperparameter-optimized model configurations are examined.

Cross-Validation and Model Evaluation

Model evaluation is performed using patient-grouped nested cross-validation:

Five outer folds
Three inner folds for hyperparameter optimization

Classification performance is assessed using:

Accuracy
F1-score
Matthews correlation coefficient (MCC)
Precision
Recall
ROC-AUC

The Technique for Order Preference by Similarity to Ideal Solution (TOPSIS) is applied to rank model configurations using multiple performance criteria.

The mature-arm workflow additionally produces feature-selection results, model comparison tables, biomarker direction information, and visualization outputs.

3. Software Environment
Python Environment

The computational environment includes:

Python 3.9.23
NumPy 2.0.2
pandas 2.3.1
scikit-learn 1.6.1
SciPy 1.13.1
imbalanced-learn 0.12.4
matplotlib 3.9.4
openpyxl 3.1.5
mrmr-selection 0.2.8

The package versions listed above define the environment used to produce the reported results.

R Environment

The R environment includes:

R 4.5.1
TCGAbiolinks 2.36.0
SummarizedExperiment 1.38.1
dplyr 1.1.4
stringr 1.5.2
tibble 3.3.0
readr 2.1.5
purrr 1.1.0

Additional R packages used in the mature-arm workflow include:

tidyr
AnnotationDbi
miRBaseVersions.db
4. Repository Structure
text
thyroid-miRNA-ML/
│
├── R/
│   └── miRNA-expression-analysis/
│       └── 01_TCGA_THCA_data_download.R
│
├── Python/
│   └── miRNA-expression-analysis/
│       └── 2209-duru-kod-Copy2-30.ipynb
│
├── mature-arm-analysis/
│   ├── 01_TCGA_THCA_mature_arm_download.R
│   └── 02_THCA_mature_arm_ML.ipynb
│
├── CITATION.cff
├── LICENSE
└── README.md
5. Reproducibility

The computational workflows are provided to support methodological transparency and reproducibility.

The mature-arm analysis was executed independently on two separate machines with separate software installations. The resulting performance tables and selected-feature sets were identical across all 180 configurations, indicating that the pipeline is deterministic under the fixed random seeds used throughout (random_state = 42).

Running the notebooks regenerates the performance tables, selected-feature sets and figures reported in the study.

miRNA Expression-Level Workflow

The corresponding R and Python scripts are available in the R/ and Python/ directories.

Mature miRNA Arm-Level Workflow

To reproduce the mature-arm analysis:

Execute 01_TCGA_THCA_mature_arm_download.R to retrieve and preprocess TCGA-THCA mature-miRNA expression data.
Generate the mature-miRNA raw-count expression matrix.
Provide the generated expression matrix as input to 02_THCA_mature_arm_ML.ipynb.
Execute the Python notebook to perform normalization, filtering, feature selection, model training, cross-validation, and TOPSIS ranking.
Examine the generated performance tables, selected features, and visualization outputs.

The computational environment is documented in the repository.

6. Data Availability

The underlying TCGA-THCA sequencing data are publicly available through the National Cancer Institute Genomic Data Commons:

https://portal.gdc.cancer.gov/

Raw sequencing data are not redistributed in this repository.

7. License

This repository is distributed under the terms specified in the LICENSE file.

8. Citation

If you use the computational workflows provided in this repository, please cite:

Vural, D., Köroğlu, S. M., & İzzetoğlu, S. (2026). TCGA-THCA miRNA Machine-Learning Analysis 
