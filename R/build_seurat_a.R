## Checkpoint: builds ONE merged Seurat object for cohort A (all 34
## specimens, author-annotated nuclei) and caches it to disk. This is the
## canonical reusable asset for cohort A going forward - QC, baseline vs
## integration comparison, and the focal-population analysis should all load
## THIS instead of re-deriving their own representation from source files.
##
## Reuses cohort_a_loaded.rds (see run_load_cohort_a.R) if present - all 34
## samples share an identical features table/gene order (verified), so
## cbind-ing their matrices needs no re-indexing.
##
## Rownames stay Ensembl gene_id (NOT symbol) - cohort A has 24 duplicate
## symbols (see QC/duplicate-check notes), so symbol can't safely be the row
## identifier. The gene_id/gene_symbol lookup table is kept in obj@misc$features
## for anything downstream that needs symbols (e.g. the ECM panel).
##
## Run: Rscript R/build_seurat_a.R
## Output: <SCRATCH_DIR>/seurat_A_raw.rds

.args <- commandArgs(trailingOnly = FALSE)
.script_path <- sub("--file=", "", .args[grep("--file=", .args)])
repo <- normalizePath(file.path(dirname(.script_path), ".."))
readRenviron(file.path(repo, ".Renviron"))
source(file.path(repo, "config.R"))
source(file.path(repo, "R", "qc_utils.R"))
library(Seurat)
library(Matrix)

cache_path <- file.path(SCRATCH_DIR, "cohort_a_loaded.rds")
if (!file.exists(cache_path)) {
  stop("cohort_a_loaded.rds not found - run R/run_load_cohort_a.R first.")
}
cat("Loading cached per-sample matrices...\n")
cohort_a <- readRDS(cache_path)

sample_manifest <- read.delim(SAMPLE_MANIFEST_PATH, stringsAsFactors = FALSE)
a_manifest <- sample_manifest[sample_manifest$cohort == "A", , drop = FALSE]

cat("Checking all samples share an identical features table (required for a direct cbind)...\n")
feat1 <- cohort_a[[1]]$features
identical_features <- all(sapply(cohort_a[-1], function(r) identical(r$features, feat1)))
if (!identical_features) {
  stop("Cohort A samples do NOT share an identical features table - cbind would silently misalign genes.")
}
cat("Confirmed identical across all", length(cohort_a), "samples.\n")

cat("Binding matrices...\n")
## Raw 10x barcodes (e.g. "AAACCCAGTTAAGAAC-1") are only unique WITHIN one
## specimen's own sequencing run - the same literal barcode string recurs
## across different specimens (34 separate 10x runs), so bare barcodes
## collide the moment matrices from different specimens are combined. Per
## DATA_DICTIONARY: identify cells by specimen + source barcode.
mat_list <- lapply(names(cohort_a), function(sid) {
  m <- cohort_a[[sid]]$matrix
  colnames(m) <- paste(sid, colnames(m), sep = "_")
  m
})
combined_matrix <- do.call(cbind, mat_list)
stopifnot(!anyDuplicated(colnames(combined_matrix)))
cat("Combined matrix:", nrow(combined_matrix), "genes x", ncol(combined_matrix), "cells\n")

cat("Building combined cell metadata...\n")
meta_list <- lapply(names(cohort_a), function(sid) {
  r <- cohort_a[[sid]]
  data.frame(
    cell_key = paste(sid, colnames(r$matrix), sep = "_"),
    sample_id = sid,
    broad_label = r$cell_metadata$Cell_type.shorter,
    detailed_label = r$cell_metadata$Cell_type.detailed,
    stringsAsFactors = FALSE
  )
})
meta <- do.call(rbind, meta_list)
meta <- merge(meta, a_manifest[, c("sample_id", "donor_id", "condition")], by = "sample_id", sort = FALSE)
meta <- meta[match(colnames(combined_matrix), meta$cell_key), ] # preserve matrix column order
rownames(meta) <- meta$cell_key
meta$cohort <- "A"

mt_ids <- mt_row_ids(feat1, id_col = "gene_id")
qc_metrics <- compute_qc_metrics(combined_matrix, mt_ids)
meta$percent_mt <- qc_metrics$percent_mt[match(meta$cell_key, qc_metrics$cell_key)]

cat("Constructing Seurat object...\n")
obj <- CreateSeuratObject(counts = combined_matrix, meta.data = meta, project = "cohort_A")
obj@misc$features <- feat1
obj@misc$cohort <- "A"

out_path <- file.path(SCRATCH_DIR, "seurat_A_raw.rds")
saveRDS(obj, out_path)
cat("\nSaved:", out_path, "\n")
cat("Dimensions:", nrow(obj), "genes x", ncol(obj), "cells\n")
print(table(obj$condition))
