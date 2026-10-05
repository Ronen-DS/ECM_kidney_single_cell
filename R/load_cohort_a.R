## Cohort A (GSE227898): raw Matrix Market archives + author annotation join.
##
## Each sample archive is <sample_id>/outs/raw_feature_bc_matrix/{matrix.mtx.gz,
## features.tsv.gz, barcodes.tsv.gz} (droplets, features x barcodes). The raw
## barcode universe is far larger than the author's analyzed cell set; we
## reconstruct that set by joining to anchor_annotations.tsv.gz on
## orig.ident (== sample_id) + individual_barcode (== raw 10x barcode).
## features.tsv.gz keeps gene ID (col 1) and symbol (col 2) both - read
## manually (not via Seurat::Read10X's symbol-only rownames) so Ensembl IDs
## survive for cross-cohort feature matching later.

library(Matrix)

#' Extracts one sample's tar.gz archive into scratch space (once; skips if
#' already extracted) and returns the raw_feature_bc_matrix directory.
extract_cohort_a_sample <- function(archive_path, sample_id) {
  dest <- file.path(SCRATCH_DIR, "cohort_A", sample_id)
  mat_dir <- file.path(dest, sample_id, "outs", "raw_feature_bc_matrix")
  if (!dir.exists(mat_dir)) {
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    utils::untar(archive_path, exdir = dest)
  }
  stopifnot(dir.exists(mat_dir))
  mat_dir
}

#' Reads one raw_feature_bc_matrix directory as a sparse dgCMatrix
#' (features x barcodes), returning the matrix plus a features data.frame
#' with gene_id/gene_symbol/feature_type (features.tsv.gz's 3 columns) -
#' kept separate from rownames so duplicate symbols never silently collide.
read_10x_raw_matrix <- function(mat_dir) {
  features <- read.delim(
    file.path(mat_dir, "features.tsv.gz"),
    header = FALSE, stringsAsFactors = FALSE,
    col.names = c("gene_id", "gene_symbol", "feature_type")
  )
  barcodes <- readLines(gzfile(file.path(mat_dir, "barcodes.tsv.gz")))
  mat <- Matrix::readMM(gzfile(file.path(mat_dir, "matrix.mtx.gz")))
  mat <- as(mat, "CsparseMatrix")
  rownames(mat) <- features$gene_id
  colnames(mat) <- barcodes
  list(matrix = mat, features = features)
}

#' Loads one cohort A specimen's counts restricted to the author's annotated
#' nuclei for that specimen, with broad/detailed labels attached.
#' `anchor_annotations` is the full table (read once by the caller, not
#' per-sample - see load_cohort_a_all()).
load_cohort_a_sample <- function(sample_id, archive_path, anchor_annotations) {
  mat_dir <- extract_cohort_a_sample(archive_path, sample_id)
  raw <- read_10x_raw_matrix(mat_dir)

  ann <- anchor_annotations[anchor_annotations$orig.ident == sample_id, , drop = FALSE]
  keep_barcodes <- intersect(colnames(raw$matrix), ann$individual_barcode)

  list(
    sample_id = sample_id,
    matrix = raw$matrix[, keep_barcodes, drop = FALSE],
    features = raw$features,
    cell_metadata = ann[match(keep_barcodes, ann$individual_barcode), , drop = FALSE],
    n_raw_barcodes = ncol(raw$matrix),
    n_annotated_in_manifest = nrow(ann),
    n_matched = length(keep_barcodes)
  )
}

#' Loads every cohort A specimen listed in the sample manifest. Returns a
#' list of per-sample results (see load_cohort_a_sample) - kept as a list
#' rather than eagerly cbind-ing, so a caller can inspect per-sample join
#' rates before deciding how to combine them.
load_cohort_a_all <- function(sample_manifest, progress = TRUE) {
  anchor_annotations <- read.delim(
    ANCHOR_ANNOTATIONS_PATH, stringsAsFactors = FALSE
  )

  a_manifest <- sample_manifest[sample_manifest$cohort == "A", , drop = FALSE]
  results <- vector("list", nrow(a_manifest))
  names(results) <- a_manifest$sample_id

  for (i in seq_len(nrow(a_manifest))) {
    sample_id <- a_manifest$sample_id[i]
    archive_path <- file.path(DATA_ROOT, a_manifest$file[i])
    if (progress) cat(sprintf("[A %d/%d] %s\n", i, nrow(a_manifest), sample_id))
    results[[sample_id]] <- load_cohort_a_sample(sample_id, archive_path, anchor_annotations)
  }
  results
}
