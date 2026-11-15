############################################################
# TCGA-THCA MATURE miRNA (3p / 5p)
# GDC Isoform Expression Quantification
#
# - Primary Tumor ve Solid Tissue Normal
# - Yalnızca "mature,MIMAT..." anotasyonlu isoformlar
# - Aynı mature MIMAT'a ait isomiR read_count'ları toplanır
# - MIMAT -> hsa-miR-xxx-3p / 5p adına çevrilir
# - Hasta + örnek tipi başına tek örnek tutulur
#   (mature-miRNA toplam library size en büyük olan)
# - Primary Tumor = 1
# - Solid Tissue Normal = 0
# - HAM COUNT çıktısı oluşturulur
# - Normalizasyon daha sonra Python'da yapılabilir
############################################################

options(stringsAsFactors = FALSE)

# ==========================================================
# Yardımcı fonksiyonlar
# ==========================================================

msg <- function(...) {
  cat("[INFO] ", paste0(...), "\n", sep = "")
}

warn <- function(...) {
  cat("[WARN] ", paste0(...), "\n", sep = "")
}

fail <- function(...) {
  stop(paste0(...), call. = FALSE)
}

install_if_missing <- function(pkgs) {
  missing <- pkgs[!(pkgs %in% rownames(installed.packages()))]
  
  if (length(missing) > 0) {
    install.packages(
      missing,
      repos = "https://cloud.r-project.org"
    )
  }
}

bioc_install_if_missing <- function(pkgs) {
  
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages(
      "BiocManager",
      repos = "https://cloud.r-project.org"
    )
  }
  
  for (p in pkgs) {
    
    if (!requireNamespace(p, quietly = TRUE)) {
      BiocManager::install(
        p,
        ask = FALSE,
        update = FALSE
      )
    }
  }
}


# ==========================================================
# Paketler
# ==========================================================

install_if_missing(
  c(
    "dplyr",
    "tidyr",
    "stringr",
    "tibble",
    "readr"
  )
)

bioc_install_if_missing(
  c(
    "TCGAbiolinks",
    "AnnotationDbi",
    "miRBaseVersions.db"
  )
)

suppressPackageStartupMessages({
  
  library(TCGAbiolinks)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(readr)
  library(AnnotationDbi)
  library(miRBaseVersions.db)
  
})


# ==========================================================
# Parametreler
# ==========================================================

project_id <- "TCGA-THCA"

out_counts <-
  "THCA_mature_miRNA_3p5p_unique_patients.csv"

out_counts_numeric <-
  "THCA_mature_miRNA_3p5p_unique_patients_NUM.csv"

out_annotation <-
  "THCA_MIMAT_to_mature_miRNA_annotation.csv"


# ==========================================================
# GDC sorgusu
# ==========================================================

msg("Isoform Expression Quantification sorgusu hazırlanıyor...")

query <- tryCatch({
  
  GDCquery(
    project = project_id,
    data.category = "Transcriptome Profiling",
    data.type = "Isoform Expression Quantification",
    experimental.strategy = "miRNA-Seq",
    sample.type = c(
      "Primary Tumor",
      "Solid Tissue Normal"
    )
  )
  
}, error = function(e) {
  
  fail(
    "GDCquery hatası: ",
    e$message
  )
  
})


# ==========================================================
# İndirme
# ==========================================================

msg("GDCdownload başlıyor...")

tryCatch({
  
  GDCdownload(
    query,
    files.per.chunk = 20
  )
  
}, error = function(e) {
  
  fail(
    "GDCdownload hatası: ",
    e$message
  )
  
})


# ==========================================================
# GDCprepare
# ==========================================================

msg("GDCprepare başlıyor...")

iso <- tryCatch({
  
  GDCprepare(
    query = query,
    summarizedExperiment = FALSE
  )
  
}, error = function(e) {
  
  fail(
    "GDCprepare hatası: ",
    e$message
  )
  
})

msg(
  "GDCprepare sınıfı: ",
  paste(
    class(iso),
    collapse = ", "
  )
)


# ==========================================================
# Gerekli sütunları kontrol et
# ==========================================================

required_cols <- c(
  "miRNA_ID",
  "read_count",
  "miRNA_region",
  "barcode"
)

missing_cols <-
  setdiff(
    required_cols,
    colnames(iso)
  )

if (length(missing_cols) > 0) {
  
  fail(
    "Eksik sütun(lar): ",
    paste(
      missing_cols,
      collapse = ", "
    )
  )
}

msg(
  "Toplam isoform satırı: ",
  nrow(iso)
)


# ==========================================================
# Yalnız mature miRNA satırlarını tut
# ==========================================================

mature <- iso |>
  
  dplyr::filter(
    stringr::str_detect(
      miRNA_region,
      "^mature,"
    )
  )

msg(
  "Mature anotasyonlu isoform satırı: ",
  nrow(mature)
)

if (nrow(mature) == 0) {
  
  fail(
    "Mature miRNA satırı bulunamadı."
  )
}


# ==========================================================
# MIMAT accession çıkar
#
# Örnek:
# mature,MIMAT0000062
# ->
# MIMAT0000062
# ==========================================================

mature <- mature |>
  
  dplyr::mutate(
    
    MIMAT =
      stringr::str_extract(
        miRNA_region,
        "MIMAT[0-9]+"
      ),
    
    read_count =
      as.numeric(
        read_count
      )
    
  ) |>
  
  dplyr::filter(
    !is.na(MIMAT)
  )


# Non-finite / NA count temizliği
mature$read_count[
  !is.finite(mature$read_count)
] <- 0

mature$read_count[
  is.na(mature$read_count)
] <- 0


# ==========================================================
# Aynı sample + MIMAT için isomiR count'larını TOPLA
#
# Bu adım kritik:
# farklı isoform koordinatları aynı mature strand'a ait olabilir.
# ==========================================================

mature_sum <- mature |>
  
  dplyr::group_by(
    barcode,
    MIMAT
  ) |>
  
  dplyr::summarise(
    
    read_count =
      sum(
        read_count,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


msg(
  "MIMAT bazında aggregate edilmiş satır sayısı: ",
  nrow(mature_sum)
)


# ==========================================================
# MIMAT -> mature miRNA adı
# miRBase v21
#
# Örn:
# MIMAT0000062 -> hsa-let-7a-5p
# ==========================================================

mimat_ids <-
  unique(
    mature_sum$MIMAT
  )

anno <- AnnotationDbi::select(
  
  miRBaseVersions.db,
  
  keys = mimat_ids,
  
  keytype = "VW-MIMAT-21.0",
  
  columns = c(
    "ACCESSION",
    "NAME"
  )
  
)


# Duplicate annotation varsa temizle
anno <- anno |>
  
  dplyr::filter(
    !is.na(ACCESSION),
    !is.na(NAME)
  ) |>
  
  dplyr::distinct(
    ACCESSION,
    .keep_all = TRUE
  )


# Annotation tablosunu kaydet
readr::write_csv(
  anno,
  out_annotation
)

msg(
  "MIMAT annotation yazıldı: ",
  out_annotation
)


# Mature isimlerini ekle
mature_sum <- mature_sum |>
  
  dplyr::left_join(
    
    anno |>
      dplyr::select(
        ACCESSION,
        NAME
      ),
    
    by = c(
      "MIMAT" = "ACCESSION"
    )
    
  )


# İsmi bulunamayanlar MIMAT olarak kalabilir
mature_sum <- mature_sum |>
  
  dplyr::mutate(
    
    mature_miRNA =
      ifelse(
        is.na(NAME),
        MIMAT,
        NAME
      )
    
  )


msg(
  "Benzersiz mature miRNA sayısı: ",
  dplyr::n_distinct(
    mature_sum$mature_miRNA
  )
)


# ==========================================================
# Aynı mature isim birden fazla precursor locus'tan
# gelebileceği için SAMPLE + MATURE NAME düzeyinde tekrar topla
#
# Bu işlem mature-arm expression matrisi oluşturur.
# ==========================================================

mature_by_name <- mature_sum |>
  
  dplyr::group_by(
    barcode,
    mature_miRNA
  ) |>
  
  dplyr::summarise(
    
    read_count =
      sum(
        read_count,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


# ==========================================================
# Sample metadata
# ==========================================================

samples <- unique(
  mature_by_name$barcode
)

meta <- tibble::tibble(
  
  barcode = samples,
  
  patient =
    substr(
      samples,
      1,
      12
    )
  
)


# TCGA sample-type kodu:
# örn:
# TCGA-XX-XXXX-01A -> 01 = Primary Tumor
# TCGA-XX-XXXX-11A -> 11 = Solid Tissue Normal

meta <- meta |>
  
  dplyr::mutate(
    
    sample_code =
      substr(
        barcode,
        14,
        15
      ),
    
    sample_type_code =
      dplyr::case_when(
        
        sample_code == "01" ~ "TP",
        
        sample_code %in%
          c(
            "11",
            "21",
            "22"
          ) ~ "NT",
        
        TRUE ~ NA_character_
      )
    
  ) |>
  
  dplyr::filter(
    sample_type_code %in%
      c(
        "TP",
        "NT"
      )
  )


msg(
  "TP/NT örnek sayısı: ",
  nrow(meta)
)


# Yalnız uygun sample'ları bırak
mature_by_name <- mature_by_name |>
  
  dplyr::filter(
    barcode %in%
      meta$barcode
  )


# ==========================================================
# Sample library size
#
# Burada mature MIMAT read count toplamı kullanılıyor.
# ==========================================================

libsizes <- mature_by_name |>
  
  dplyr::group_by(
    barcode
  ) |>
  
  dplyr::summarise(
    
    libsize =
      sum(
        read_count,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


meta <- meta |>
  
  dplyr::left_join(
    libsizes,
    by = "barcode"
  )


# ==========================================================
# Hasta + örnek tipi başına tek örnek
# En yüksek mature-miRNA library size
# ==========================================================

meta_unique <- meta |>
  
  dplyr::group_by(
    patient,
    sample_type_code
  ) |>
  
  dplyr::slice_max(
    order_by = libsize,
    n = 1,
    with_ties = FALSE
  ) |>
  
  dplyr::ungroup()


msg(
  "Duplicate çıkarıldıktan sonra örnek sayısı: ",
  nrow(meta_unique)
)


# ==========================================================
# Label
# ==========================================================

meta_unique <- meta_unique |>
  
  dplyr::mutate(
    
    label =
      dplyr::case_when(
        
        sample_type_code == "TP" ~ 1L,
        
        sample_type_code == "NT" ~ 0L,
        
        TRUE ~ NA_integer_
      )
    
  )


lab_tab <-
  table(
    meta_unique$label
  )

msg(
  "Label dağılımı (0=Normal, 1=Tumor): ",
  paste(
    names(lab_tab),
    lab_tab,
    collapse = " | "
  )
)


# ==========================================================
# Final sample set
# ==========================================================

keep_samples <-
  meta_unique$barcode

mature_final <- mature_by_name |>
  
  dplyr::filter(
    barcode %in%
      keep_samples
  )


# ==========================================================
# LONG -> WIDE
#
# Sample
# label
# hsa-let-7a-5p
# hsa-let-7a-3p
# hsa-miR-...
# ==========================================================

count_matrix <- mature_final |>
  
  dplyr::select(
    barcode,
    mature_miRNA,
    read_count
  ) |>
  
  tidyr::pivot_wider(
    
    names_from =
      mature_miRNA,
    
    values_from =
      read_count,
    
    values_fill = 0
    
  )


# ==========================================================
# Label ekle
# ==========================================================

ml_df <- count_matrix |>
  
  dplyr::left_join(
    
    meta_unique |>
      dplyr::select(
        barcode,
        label
      ),
    
    by = "barcode"
    
  ) |>
  
  dplyr::rename(
    Sample = barcode
  ) |>
  
  dplyr::relocate(
    Sample,
    label
  )


# ==========================================================
# Kontroller
# ==========================================================

if (nrow(ml_df) == 0) {
  fail(
    "Final mature-miRNA matrisi boş."
  )
}

if (any(is.na(ml_df$label))) {
  fail(
    "Label sütununda NA var."
  )
}


feature_cols <-
  setdiff(
    names(ml_df),
    c(
      "Sample",
      "label"
    )
  )


ml_df[
  feature_cols
] <- lapply(
  
  ml_df[
    feature_cols
  ],
  
  function(x) {
    
    x <- as.numeric(x)
    
    x[
      is.na(x) |
        !is.finite(x)
    ] <- 0
    
    x
  }
  
)


# ==========================================================
# Sample + label + mature-miRNA HAM COUNTS
# ==========================================================

readr::write_csv(
  ml_df,
  out_counts
)

msg(
  "Yazıldı: ",
  out_counts
)


# ==========================================================
# Numeric-only çıktı
# ==========================================================

num_df <- ml_df |>
  
  dplyr::select(
    -Sample
  ) |>
  
  dplyr::relocate(
    label
  )


readr::write_csv(
  num_df,
  out_counts_numeric
)

msg(
  "Yazıldı: ",
  out_counts_numeric
)


# ==========================================================
# Son özet
# ==========================================================

msg(
  "Final örnek sayısı: ",
  nrow(ml_df)
)

msg(
  "Final mature-miRNA feature sayısı: ",
  length(feature_cols)
)

msg(
  "Tamamlandı."
)

msg(
  "Çıktılar:"
)

msg(
  " - ",
  out_counts
)

msg(
  " - ",
  out_counts_numeric
)

msg(
  " - ",
  out_annotation
)
