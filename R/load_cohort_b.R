## Cohort B (GSE178481): compressed count CSVs (genes x cells, gene symbols
## as row names) inside one tar archive, joined to cohort_b_annotations.csv
## on the cell key (specimen name + "_" + barcode - already how this
## archive's CSV column headers are formatted, confirmed by inspection; no
## extra prefixing needed for these files specifically).
##
## B has no Ensembl IDs at all (symbols only) - duplicate symbols are
## reported here, not silently collapsed; a cross-cohort harmonization step
## decides how to handle them.

#' Extracts one sample's count.csv.gz member from the shared tar archive
#' into scratch space (once; skips if already extracted).
extract_cohort_b_member <- function(archive_member) {
  dest_dir <- file.path(SCRATCH_DIR, "cohort_B")
  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
  csv_path <- file.path(dest_dir, archive_member)
  if (!file.exists(csv_path)) {
    tar_path <- file.path(COHORT_B_DIR, "GSE178481_RAW.tar")
    utils::untar(tar_path, exdir = dest_dir, files = archive_member)
  }
  stopifnot(file.exists(csv_path))
  csv_path
}

#' Loads one cohort B specimen: genes (symbol) x cells sparse matrix, plus
#' cell metadata from cohort_b_annotations.csv (ano.l1/ano.l2 labels and the
#' source QC fields) for whichever of this specimen's columns are annotated.
load_cohort_b_sample <- function(sample_id, archive_member, annotations) {
  csv_path <- extract_cohort_b_member(archive_member)

  df <- read.csv(csv_path, row.names = 1, check.names = FALSE)
  mat <- as(as.matrix(df), "CsparseMatrix")
  # read.csv mangles a leading digit/symbol in column names (e.g. prepends
  # "X"); count.csv.gz columns are already valid, well-formed cell keys
  # (specimen_barcode), so restore them exactly as they appear in the file.
  colnames(mat) <- colnames(df)

  ann <- annotations[annotations$cell_key %in% colnames(mat), , drop = FALSE]
  keep <- intersect(colnames(mat), ann$cell_key)

  dup_symbols <- rownames(mat)[duplicated(rownames(mat))]

  list(
    sample_id = sample_id,
    matrix = mat[, keep, drop = FALSE],
    cell_metadata = ann[match(keep, ann$cell_key), , drop = FALSE],
    n_raw_cells = ncol(mat),
    n_annotated_in_file = sum(annotations$cell_key %in% colnames(mat)),
    n_matched = length(keep),
    duplicate_gene_symbols = unique(dup_symbols)
  )
}

#' Loads every cohort B specimen listed in the sample manifest.
load_cohort_b_all <- function(sample_manifest, progress = TRUE) {
  annotations <- read.csv(COHORT_B_ANNOTATIONS_PATH, stringsAsFactors = FALSE)
  colnames(annotations)[colnames(annotations) == "X"] <- "cell_key"

  b_manifest <- sample_manifest[sample_manifest$cohort == "B", , drop = FALSE]
  results <- vector("list", nrow(b_manifest))
  names(results) <- b_manifest$sample_id

  for (i in seq_len(nrow(b_manifest))) {
    sample_id <- b_manifest$sample_id[i]
    if (progress) cat(sprintf("[B %d/%d] %s\n", i, nrow(b_manifest), sample_id))
    results[[sample_id]] <- load_cohort_b_sample(
      sample_id, b_manifest$archive_member[i], annotations
    )
  }
  results
}
