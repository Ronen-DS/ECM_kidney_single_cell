## Collapses Ensembl-ID duplicate rows (multiple ENSG IDs mapping to the same
## gene symbol - a GENCODE reference quirk affecting cohorts A and C, see the
## duplicate-ID audit) down to one row per gene symbol, by summing counts.
## Cohort B's matrix is already symbol-indexed with zero duplicates, so it
## passes through unchanged (still re-saved under the same *_symbols.rds
## name, so all 3 cohorts have the same downstream entry point).
##
## Input:  seurat_{A,B,C}_filtered.rds (post-QC-filter checkpoints)
## Output: seurat_{A,B,C}_symbols.rds  (same cells, symbol-indexed features)
##         results/tables/duplicate_symbol_summary.tsv
##         results/tables/duplicate_symbol_summary_by_cohort.tsv
##
## Written so this can later be pasted into notebook cells as-is: no
## commandArgs-only path logic - falls back to the working directory as repo
## root when not invoked via `Rscript path\to\file.R`.
##
## Run: Rscript R/collapse_duplicate_symbols.R

repo <- local({
  args <- commandArgs(trailingOnly = FALSE)
  script_arg <- args[grep("^--file=", args)]
  if (length(script_arg) == 1) {
    normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."))
  } else {
    ## Not run via Rscript (e.g. a notebook cell or interactive session) -
    ## assume the working directory is the repo root.
    normalizePath(".")
  }
})
readRenviron(file.path(repo, ".Renviron"))
source(file.path(repo, "config.R"))
source(file.path(repo, "R", "symbol_utils.R"))
library(Seurat)
library(Matrix)

## Drops meta.data's own nCount_RNA/nFeature_RNA (computed on the pre-collapse,
## Ensembl-ID-level matrix) so CreateSeuratObject recomputes them fresh from
## the symbol-collapsed matrix instead of silently reusing stale values.
## nCount_RNA is unaffected by collapsing (row sums are conserved by
## summing), but nFeature_RNA can drop slightly wherever a cell had nonzero
## counts in more than one now-merged duplicate row.
drop_stale_qc_cols <- function(meta) {
  meta$nCount_RNA <- NULL
  meta$nFeature_RNA <- NULL
  meta
}

## ---- Cohort A ---------------------------------------------------------------
cat("Cohort A...\n")
obj_a <- readRDS(file.path(SCRATCH_DIR, "seurat_A_filtered.rds"))
features_a <- obj_a@misc$features
counts_a <- GetAssayData(obj_a, assay = "RNA", layer = "counts")
new_counts_a <- collapse_symbol_duplicates(counts_a, features_a)
obj_a_symbols <- CreateSeuratObject(counts = new_counts_a, meta.data = drop_stale_qc_cols(obj_a@meta.data), project = "cohort_A")
obj_a_symbols@misc$cohort <- "A"
obj_a_symbols@misc$features_pre_collapse <- features_a
saveRDS(obj_a_symbols, file.path(SCRATCH_DIR, "seurat_A_symbols.rds"))
summary_a <- duplicate_symbol_table(features_a, "A")
cat(sprintf("  %d Ensembl IDs -> %d symbols (%d duplicated symbols, %d extra rows summed)\n",
            nrow(features_a), nrow(obj_a_symbols), nrow(summary_a),
            sum(summary_a$n_duplicates) - nrow(summary_a)))

## ---- Cohort B ---------------------------------------------------------------
cat("Cohort B...\n")
obj_b <- readRDS(file.path(SCRATCH_DIR, "seurat_B_filtered.rds"))
stopifnot(!anyDuplicated(rownames(obj_b))) # already symbol-indexed, verified zero duplicates
saveRDS(obj_b, file.path(SCRATCH_DIR, "seurat_B_symbols.rds"))
summary_b <- duplicate_symbol_table(data.frame(gene_id = character(0), gene_symbol = character(0)), "B")
cat(sprintf("  %d genes, already symbol-indexed, 0 duplicates - passthrough\n", nrow(obj_b)))

## ---- Cohort C ---------------------------------------------------------------
cat("Cohort C...\n")
obj_c <- readRDS(file.path(SCRATCH_DIR, "seurat_C_filtered.rds"))
features_c <- obj_c@misc$features
counts_c <- GetAssayData(obj_c, assay = "RNA", layer = "counts")
new_counts_c <- collapse_symbol_duplicates(counts_c, features_c)
obj_c_symbols <- CreateSeuratObject(counts = new_counts_c, meta.data = drop_stale_qc_cols(obj_c@meta.data), project = "cohort_C")
obj_c_symbols@misc$cohort <- "C"
obj_c_symbols@misc$features_pre_collapse <- features_c
saveRDS(obj_c_symbols, file.path(SCRATCH_DIR, "seurat_C_symbols.rds"))
summary_c <- duplicate_symbol_table(features_c, "C")
cat(sprintf("  %d Ensembl IDs -> %d symbols (%d duplicated symbols, %d extra rows summed)\n",
            nrow(features_c), nrow(obj_c_symbols), nrow(summary_c),
            sum(summary_c$n_duplicates) - nrow(summary_c)))

## ---- Summary tables -----------------------------------------------------------
summary_all <- rbind(summary_a, summary_b, summary_c)
tables_dir <- file.path(repo, "results", "tables")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
write.table(summary_all, file.path(tables_dir, "duplicate_symbol_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

per_cohort <- data.frame(
  cohort = c("A", "B", "C"),
  n_duplicated_symbols = c(nrow(summary_a), nrow(summary_b), nrow(summary_c)),
  n_extra_rows_summed = c(sum(summary_a$n_duplicates) - nrow(summary_a),
                           sum(summary_b$n_duplicates) - nrow(summary_b),
                           sum(summary_c$n_duplicates) - nrow(summary_c))
)
write.table(per_cohort, file.path(tables_dir, "duplicate_symbol_summary_by_cohort.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nPer-cohort duplicate-symbol counts:\n")
print(per_cohort, row.names = FALSE)
cat(sprintf("\nWrote %d rows to results/tables/duplicate_symbol_summary.tsv\n", nrow(summary_all)))
