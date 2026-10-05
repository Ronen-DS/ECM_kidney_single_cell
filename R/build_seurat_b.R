## Checkpoint: builds ONE merged Seurat object for cohort B (all 26
## specimens) and caches it to disk - the canonical reusable asset for
## cohort B going forward (QC, integration, focal-population analysis all
## load THIS instead of re-parsing the 26 raw count.csv.gz files, which is
## the expensive step here: ~35 min, dominated by a few large dense CSV
## reads for the biggest specimens).
##
## Unlike cohort A/C, B's raw column headers already embed the specimen name
## (e.g. "RCC-PR6-Normal_AAACGGGAGACTAAGT-1" - a property of the source
## files themselves, not something added here), so no extra barcode
## prefixing is needed before combining - just a safety check that it's
## actually true.
##
## Rownames stay gene SYMBOL (B has no Ensembl ID at all). B's supplied
## author QC fields (TotalUMI, gene.num, MT.ratio, doublest.score) and
## broad/detailed labels (ano.l1/ano.l2) are kept as metadata columns
## alongside our own recomputed percent_mt, so qc_cohort_b.R's consistency
## check can run directly off this object without re-deriving anything.
##
## Run: Rscript R/build_seurat_b.R
## Output: <SCRATCH_DIR>/seurat_B_raw.rds

.args <- commandArgs(trailingOnly = FALSE)
.script_path <- sub("--file=", "", .args[grep("--file=", .args)])
repo <- normalizePath(file.path(dirname(.script_path), ".."))
readRenviron(file.path(repo, ".Renviron"))
source(file.path(repo, "config.R"))
source(file.path(repo, "R", "load_cohort_b.R"))
source(file.path(repo, "R", "qc_utils.R"))
library(Seurat)
library(Matrix)

sample_manifest <- read.delim(SAMPLE_MANIFEST_PATH, stringsAsFactors = FALSE)
b_manifest <- sample_manifest[sample_manifest$cohort == "B", , drop = FALSE]

cat("Loading all 26 cohort B specimens (this is the ~35 min step)...\n")
cohort_b <- load_cohort_b_all(sample_manifest)

feat1 <- rownames(cohort_b[[1]]$matrix)
identical_genes <- all(sapply(cohort_b[-1], function(r) identical(rownames(r$matrix), feat1)))
if (!identical_genes) {
  stop("Cohort B samples do NOT share an identical gene list/order - cbind would silently misalign genes.")
}
cat("Confirmed identical gene list across all", length(cohort_b), "samples.\n")

mat_list <- lapply(cohort_b, function(r) r$matrix)
combined_matrix <- do.call(cbind, mat_list)
stopifnot(!anyDuplicated(colnames(combined_matrix)))
cat("Combined matrix:", nrow(combined_matrix), "genes x", ncol(combined_matrix), "cells\n")

cat("Building combined cell metadata...\n")
meta_list <- lapply(names(cohort_b), function(sid) {
  r <- cohort_b[[sid]]
  df <- r$cell_metadata[, c("cell_key", "ano.l1", "ano.l2", "TotalUMI", "doublest.score", "MT.ratio", "gene.num")]
  df$sample_id <- sid
  df
})
meta <- do.call(rbind, meta_list)
meta <- merge(meta, b_manifest[, c("sample_id", "donor_id", "condition")], by = "sample_id", sort = FALSE)
meta <- meta[match(colnames(combined_matrix), meta$cell_key), ]
rownames(meta) <- meta$cell_key
meta$cohort <- "B"

mt_ids <- grep("^MT-", rownames(combined_matrix), value = TRUE)
qc_metrics <- compute_qc_metrics(combined_matrix, mt_ids)
meta$percent_mt <- qc_metrics$percent_mt[match(meta$cell_key, qc_metrics$cell_key)]
meta$MT.ratio_pct <- meta$MT.ratio * 100

cat("Constructing Seurat object...\n")
obj <- CreateSeuratObject(counts = combined_matrix, meta.data = meta, project = "cohort_B")
obj@misc$cohort <- "B"

out_path <- file.path(SCRATCH_DIR, "seurat_B_raw.rds")
saveRDS(obj, out_path)
cat("\nSaved:", out_path, "\n")
cat("Dimensions:", nrow(obj), "genes x", ncol(obj), "cells\n")
print(table(obj$condition))
