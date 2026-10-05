## Checkpoint: builds ONE merged Seurat object for cohort C (all 13
## specimens). Fast (H5 reads, no author annotation join needed), so this
## just loads fresh each time rather than needing its own separate raw cache.
##
## Rownames stay Ensembl gene_id (cohort C has 34 duplicate symbols - see
## QC/duplicate-check notes). Features lookup kept in obj@misc$features.
##
## Run: Rscript R/build_seurat_c.R
## Output: <SCRATCH_DIR>/seurat_C_raw.rds

.args <- commandArgs(trailingOnly = FALSE)
.script_path <- sub("--file=", "", .args[grep("--file=", .args)])
repo <- normalizePath(file.path(dirname(.script_path), ".."))
readRenviron(file.path(repo, ".Renviron"))
source(file.path(repo, "config.R"))
source(file.path(repo, "R", "load_cohort_c.R"))
source(file.path(repo, "R", "qc_utils.R"))
library(Seurat)
library(Matrix)

sample_manifest <- read.delim(SAMPLE_MANIFEST_PATH, stringsAsFactors = FALSE)
c_manifest <- sample_manifest[sample_manifest$cohort == "C", , drop = FALSE]

cat("Loading all 13 cohort C specimens...\n")
cohort_c <- load_cohort_c_all(sample_manifest)

feat1 <- cohort_c[[1]]$features
identical_features <- all(sapply(cohort_c[-1], function(r) identical(r$features, feat1)))
if (!identical_features) {
  stop("Cohort C samples do NOT share an identical features table - cbind would silently misalign genes.")
}
cat("Confirmed identical features across all", length(cohort_c), "samples.\n")

## Raw barcodes are only unique WITHIN one specimen's own run - see
## build_seurat_a.R for why bare barcodes must be prefixed before combining
## matrices from different specimens.
mat_list <- lapply(names(cohort_c), function(sid) {
  m <- cohort_c[[sid]]$matrix
  colnames(m) <- paste(sid, colnames(m), sep = "_")
  m
})
combined_matrix <- do.call(cbind, mat_list)
stopifnot(!anyDuplicated(colnames(combined_matrix)))
cat("Combined matrix:", nrow(combined_matrix), "genes x", ncol(combined_matrix), "cells\n")

meta_list <- lapply(names(cohort_c), function(sid) {
  r <- cohort_c[[sid]]
  data.frame(cell_key = paste(sid, colnames(r$matrix), sep = "_"), sample_id = sid, stringsAsFactors = FALSE)
})
meta <- do.call(rbind, meta_list)
meta <- merge(meta, c_manifest[, c("sample_id", "donor_id", "condition")], by = "sample_id", sort = FALSE)
meta <- meta[match(colnames(combined_matrix), meta$cell_key), ]
rownames(meta) <- meta$cell_key
meta$cohort <- "C"

mt_ids <- mt_row_ids(feat1, id_col = "gene_id")
qc_metrics <- compute_qc_metrics(combined_matrix, mt_ids)
meta$percent_mt <- qc_metrics$percent_mt[match(meta$cell_key, qc_metrics$cell_key)]

cat("Constructing Seurat object...\n")
obj <- CreateSeuratObject(counts = combined_matrix, meta.data = meta, project = "cohort_C")
obj@misc$features <- feat1
obj@misc$cohort <- "C"

out_path <- file.path(SCRATCH_DIR, "seurat_C_raw.rds")
saveRDS(obj, out_path)
cat("\nSaved:", out_path, "\n")
cat("Dimensions:", nrow(obj), "genes x", ncol(obj), "cells\n")
print(table(obj$condition))
