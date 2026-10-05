## Step 1 of the Endothelial ECM pseudobulk Tumor-vs-Adjacent comparison:
## builds per-cohort pseudobulk (gene x sample) raw-count matrices by summing
## RNA-assay COUNTS (never the integrated assay, never an embedding) over
## Endothelial cells, grouped by donor_id x condition.
##
## Gene identifier space: seurat_{A,B,C}_symbols.rds already hold raw counts
## on a shared gene-SYMBOL row space (see R/symbol_utils.R::collapse_symbol_
## duplicates - cohort A/C's Ensembl-ID duplicate rows per symbol were summed
## there; cohort B has no Ensembl ID at all and is natively symbol-indexed).
## That IS this project's existing Ensembl-ID gene-matching step, so no
## further ID reconciliation is done here - these checkpoints are used
## directly as each cohort's full, post-QC-filter raw-count gene universe.
##
## Donor grouping: donor_id + condition (from merged_ABC_annotated.rds's own
## metadata, itself derived from metadata/sample_manifest.tsv) already pools
## repeated specimens of the same donor+condition - e.g. cohort C's Patient
## SS_2014 Adjacent (Cortex + Medulla, two separate GSM specimens per
## metadata/original/cohort_C_samples.tsv's tissue column) sum into ONE
## pseudobulk sample, with no extra join needed (see
## R/ecm_population_choice_coverage.R, which verified and documented this
## same donor+condition pooling).
##
## Samples with < 20 Endothelial cells are dropped (per task spec).
##
## Input:  merged_ABC_annotated.rds (metadata only: cohort, donor_id,
##         condition, cell_annotation), seurat_{A,B,C}_symbols.rds (raw counts)
## Output: <SCRATCH_DIR>/pseudobulk_endothelial.rds
##           list(counts = list(A=genes x samples matrix, B=..., C=...),
##                sample_table = data.frame(cohort, donor_id, condition,
##                                           sample_name, n_cells, library_size))
##         results/ecm_endothelial/pseudobulk_sample_table.tsv
##
## Run: Rscript R/ecm_endothelial_pseudobulk_build.R

repo <- local({
  args <- commandArgs(trailingOnly = FALSE)
  script_arg <- args[grep("^--file=", args)]
  if (length(script_arg) == 1) {
    normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."))
  } else {
    normalizePath(".")
  }
})
readRenviron(file.path(repo, ".Renviron"))
source(file.path(repo, "config.R"))
library(Seurat)
library(Matrix)

MIN_CELLS <- 20
out_dir <- file.path(repo, "results", "ecm_endothelial")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading merged_ABC_annotated.rds for metadata only...\n")
obj_meta <- readRDS(file.path(SCRATCH_DIR, "merged_ABC_annotated.rds"))
full_meta <- obj_meta@meta.data
full_meta$cell_key <- rownames(full_meta)
rm(obj_meta); gc(verbose = FALSE)
stopifnot(all(c("cohort", "donor_id", "condition", "cell_annotation") %in% colnames(full_meta)))

endo_meta <- full_meta[full_meta$cell_annotation == "Endothelial", ]
cat(sprintf("Endothelial cells total: %d (A=%d, B=%d, C=%d)\n", nrow(endo_meta),
            sum(endo_meta$cohort == "A"), sum(endo_meta$cohort == "B"), sum(endo_meta$cohort == "C")))

pb_counts <- list()
sample_rows <- list()

for (co in c("A", "B", "C")) {
  cat(sprintf("\nCohort %s: loading seurat_%s_symbols.rds (raw counts)...\n", co, co))
  obj <- readRDS(file.path(SCRATCH_DIR, sprintf("seurat_%s_symbols.rds", co)))
  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  rm(obj); gc(verbose = FALSE)

  sub <- endo_meta[endo_meta$cohort == co, ]
  missing <- setdiff(sub$cell_key, colnames(counts))
  if (length(missing) > 0) {
    stop(sprintf("Cohort %s: %d Endothelial cell_key(s) from merged_ABC_annotated.rds not found in seurat_%s_symbols.rds colnames - cell identifiers diverged, aborting.",
                 co, length(missing), co))
  }
  cat(sprintf("  Confirmed all %d cohort-%s Endothelial cell_keys match seurat_%s_symbols.rds colnames.\n", nrow(sub), co, co))

  counts_sub <- counts[, sub$cell_key, drop = FALSE]
  rm(counts); gc(verbose = FALSE)

  group <- paste(sub$donor_id, sub$condition, sep = "\r")
  group_levels <- unique(group)
  indicator <- sparseMatrix(i = seq_along(group), j = match(group, group_levels), x = 1,
                             dims = c(length(group), length(group_levels)))
  pb <- as(counts_sub %*% indicator, "CsparseMatrix")
  colnames(pb) <- group_levels
  rm(counts_sub); gc(verbose = FALSE)

  n_cells_per_group <- as.integer(table(factor(group, levels = group_levels)))
  lib_size <- Matrix::colSums(pb)

  keep <- n_cells_per_group >= MIN_CELLS
  cat(sprintf("  %d/%d donor x condition groups have >= %d Endothelial cells\n",
              sum(keep), length(group_levels), MIN_CELLS))

  pb_counts[[co]] <- pb[, keep, drop = FALSE]

  split_names <- strsplit(group_levels[keep], "\r", fixed = TRUE)
  sample_rows[[co]] <- data.frame(
    cohort = co,
    donor_id = vapply(split_names, `[`, character(1), 1),
    condition = vapply(split_names, `[`, character(1), 2),
    n_cells = n_cells_per_group[keep],
    library_size = lib_size[keep],
    stringsAsFactors = FALSE
  )
  rm(pb); gc(verbose = FALSE)
}

sample_table <- do.call(rbind, sample_rows)
sample_table$sample_name <- paste(sample_table$cohort, sample_table$donor_id, sample_table$condition, sep = "_")
rownames(sample_table) <- NULL
for (co in names(pb_counts)) {
  stopifnot(identical(colnames(pb_counts[[co]]),
                       paste(sample_table$donor_id[sample_table$cohort == co],
                             sample_table$condition[sample_table$cohort == co], sep = "\r")))
  colnames(pb_counts[[co]]) <- sample_table$sample_name[sample_table$cohort == co]
}

sample_table <- sample_table[order(sample_table$cohort, sample_table$donor_id, sample_table$condition), ]
write.table(sample_table, file.path(out_dir, "pseudobulk_sample_table.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== Pseudobulk sample table (Endothelial, >= 20 cells) ===\n")
print(sample_table, row.names = FALSE)
cat(sprintf("\nTotal samples: %d (A=%d, B=%d, C=%d)\n", nrow(sample_table),
            sum(sample_table$cohort == "A"), sum(sample_table$cohort == "B"), sum(sample_table$cohort == "C")))
cat("Gene universe sizes (rows) per cohort:\n")
for (co in names(pb_counts)) cat(sprintf("  %s: %d genes\n", co, nrow(pb_counts[[co]])))

ecm_genes <- c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1",
               "COL6A2", "FN1", "SPARC", "DCN", "LUM")
cat("\nECM panel gene presence in each cohort's pseudobulk gene universe:\n")
for (co in names(pb_counts)) {
  present <- ecm_genes %in% rownames(pb_counts[[co]])
  cat(sprintf("  %s: %s\n", co, paste(sprintf("%s=%s", ecm_genes, present), collapse = ", ")))
}

saveRDS(list(counts = pb_counts, sample_table = sample_table),
        file.path(SCRATCH_DIR, "pseudobulk_endothelial.rds"))
cat("\nSaved:", file.path(SCRATCH_DIR, "pseudobulk_endothelial.rds"), "\n")
cat("Done.\n")
