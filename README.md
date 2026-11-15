# TCGA-THCA: mature miRNA arm-level (3p/5p) analysis

This folder contains the mature-arm extension of the TCGA-THCA miRNA machine-learning workflow. The original analysis files can remain unchanged in the parent repository.

## Files and order

1. `01_TCGA_THCA_mature_arm_download.R` — retrieves TCGA-THCA miRNA-Seq **Isoform Expression Quantification** data for primary tumors and solid-tissue normals; keeps `mature,MIMAT...` isoforms; sums isomiR read counts by sample/MIMAT; maps MIMAT accessions to mature miRNA names using `miRBaseVersions.db` (miRBase v21); aggregates by mature miRNA name; selects one sample per patient and sample type (highest mature-miRNA library size); exports raw counts.
2. `02_THCA_mature_arm_ML.ipynb` — reads the raw counts CSV, computes RPM and `log2(RPM + 1)`, filters low-expression mature miRNAs (RPM > 1 in at least 20% of samples), and evaluates classifiers with ANOVA/MI/mRMR feature selection (`k=30`), class-balance scenarios and patient-grouped cross-validation. It also writes feature-selection and TOPSIS summaries.

## Run

1. Run the R script in an R environment with **TCGAbiolinks**, **AnnotationDbi**, **miRBaseVersions.db**, **dplyr**, **tidyr**, **stringr**, **tibble**, and **readr** available. Some packages are installed automatically if missing.
2. Keep the resulting `THCA_mature_miRNA_3p5p_unique_patients.csv` in **this folder**, or adjust `file_path` at the start of the notebook to point to its location.
3. Start Jupyter with this folder as the working directory; then run `02_THCA_mature_arm_ML.ipynb`.

Python packages used include `pandas`, `numpy`, `matplotlib`, `openpyxl`, `scikit-learn`, `imbalanced-learn`, and `mrmr-selection` (imported as `mrmr`). The notebook defines a 5-fold outer and 3-fold inner `StratifiedGroupKFold` for optimized models. Its scenarios are unbalanced, random undersampling, and SMOTE; its selectors are ANOVA, mutual information, and mRMR.

## Data outputs

The R script writes:

- `THCA_mature_miRNA_3p5p_unique_patients.csv` — sample ID, binary label and raw counts.
- `THCA_mature_miRNA_3p5p_unique_patients_NUM.csv` — numeric-only matrix (no sample IDs).
- `THCA_MIMAT_to_mature_miRNA_annotation.csv` — MIMAT/name annotation table.

`label=1` indicates primary tumor; `label=0` indicates normal tissue.

The Python notebook creates Excel results, feature-selection summaries, and PNG plots in its working directory. Input GDC data and generated results are not included in this package.

## Reproducibility note

The notebook has been checked for Python syntax, but this package has **not** been fully executed or scientifically validated against the revised manuscript. Re-run the workflow to confirm exact sample counts, retained features and result metrics. Preserve software/package versions and the GDC data download date for publication-grade reproducibility.
