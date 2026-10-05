## Mesenchymal-contamination check for the Endothelial label: flags any
## Endothelial-labeled cell carrying >=1 raw count of a canonical pericyte/
## fibroblast/smooth-muscle marker NOT in the ECM panel, then rebuilds the
## Endothelial pseudobulk from the UNFLAGGED cells only, to test whether the
## earlier Tumor-vs-Adjacent ECM signal (strong in B, null in C, opposite
## direction in C's paired donors) survives once likely mesenchymal doublet/
## contaminant cells are removed from the Endothelial label.
##
## Markers: RGS5, PDGFRB, ACTA2, PDGFRA, MYH11, NOTCH3 - confirmed present by
## exact symbol match in ALL THREE cohorts' own seurat_{A,B,C}_symbols.rds
## gene universe (same checkpoint the pseudobulk itself is built from), via
## R/check scripts run before this one. None of these were "absent-cohort"
## cases like PECAM1 in B, so no cohort needs special-casing here. DCN, LUM
## and all collagens are deliberately excluded - they are the outcome being
## tested, not a contamination marker.
##
## Exclusion rule, not inclusion: a cell is flagged if it has >=1 count of
## ANY marker, applied identically to all 3 cohorts. An inclusion rule
## ("must express an endothelial marker") was rejected per the task spec,
## since cohort C's heavier marker dropout would disproportionately remove
## real C endothelial cells under that rule, confounding the comparison.
##
## Input:  merged_ABC_annotated.rds (metadata), seurat_{A,B,C}_symbols.rds (raw counts)
## Output: results/ecm_endothelial/mesenchymal_flag_summary.tsv (per cohort x condition)
##         results/ecm_endothelial/mesenchymal_flag_cells.tsv (per cell)
##         results/ecm_endothelial/pseudobulk_clean_sample_table.tsv
##         <SCRATCH_DIR>/pseudobulk_endothelial_clean.rds
##
## Run: Rscript R/ecm_endothelial_mesenchymal_flag.R

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
mes_markers <- c("RGS5", "PDGFRB", "ACTA2", "PDGFRA", "MYH11", "NOTCH3")
out_dir <- file.path(repo, "results", "ecm_endothelial")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading merged_ABC_annotated.rds for metadata only...\n")
obj_meta <- readRDS(file.path(SCRATCH_DIR, "merged_ABC_annotated.rds"))
full_meta <- obj_meta@meta.data
full_meta$cell_key <- rownames(full_meta)
rm(obj_meta); gc(verbose = FALSE)
endo_meta <- full_meta[full_meta$cell_annotation == "Endothelial", ]
cat(sprintf("Endothelial cells total: %d (A=%d, B=%d, C=%d)\n", nrow(endo_meta),
            sum(endo_meta$cohort == "A"), sum(endo_meta$cohort == "B"), sum(endo_meta$cohort == "C")))

flag_rows <- list()
pb_counts_clean <- list()
sample_rows_clean <- list()

for (co in c("A", "B", "C")) {
  cat(sprintf("\nCohort %s: loading seurat_%s_symbols.rds (raw counts)...\n", co, co))
  obj <- readRDS(file.path(SCRATCH_DIR, sprintf("seurat_%s_symbols.rds", co)))
  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  rm(obj); gc(verbose = FALSE)

  present_markers <- intersect(mes_markers, rownames(counts))
  cat(sprintf("  Mesenchymal markers present: %s\n", paste(present_markers, collapse = ", ")))
  stopifnot(length(present_markers) == length(mes_markers))  # all 6 confirmed present in all 3 cohorts beforehand

  sub <- endo_meta[endo_meta$cohort == co, ]
  stopifnot(all(sub$cell_key %in% colnames(counts)))
  counts_endo <- counts[, sub$cell_key, drop = FALSE]
  rm(counts); gc(verbose = FALSE)

  marker_counts <- counts_endo[present_markers, , drop = FALSE]
  flagged <- Matrix::colSums(marker_counts) >= 1
  cat(sprintf("  Flagged (>=1 mesenchymal marker count): %d/%d (%.1f%%)\n",
              sum(flagged), length(flagged), 100 * mean(flagged)))

  flag_rows[[co]] <- data.frame(cell_key = sub$cell_key, cohort = co, donor_id = sub$donor_id,
                                 condition = sub$condition, flagged = flagged, stringsAsFactors = FALSE)

  ## --- Clean pseudobulk: unflagged cells only, same donor x condition summation ---
  clean_cells <- sub$cell_key[!flagged]
  clean_meta <- sub[!flagged, ]
  counts_clean <- counts_endo[, clean_cells, drop = FALSE]
  rm(counts_endo); gc(verbose = FALSE)

  group <- paste(clean_meta$donor_id, clean_meta$condition, sep = "\r")
  group_levels <- unique(group)
  indicator <- sparseMatrix(i = seq_along(group), j = match(group, group_levels), x = 1,
                             dims = c(length(group), length(group_levels)))
  pb <- as(counts_clean %*% indicator, "CsparseMatrix")
  colnames(pb) <- group_levels
  rm(counts_clean); gc(verbose = FALSE)

  n_cells_per_group <- as.integer(table(factor(group, levels = group_levels)))
  lib_size <- Matrix::colSums(pb)
  keep <- n_cells_per_group >= MIN_CELLS
  cat(sprintf("  Clean pseudobulk: %d/%d donor x condition groups have >= %d unflagged Endothelial cells\n",
              sum(keep), length(group_levels), MIN_CELLS))

  pb_counts_clean[[co]] <- pb[, keep, drop = FALSE]
  split_names <- strsplit(group_levels[keep], "\r", fixed = TRUE)
  sample_rows_clean[[co]] <- data.frame(
    cohort = co, donor_id = vapply(split_names, `[`, character(1), 1),
    condition = vapply(split_names, `[`, character(1), 2),
    n_cells = n_cells_per_group[keep], library_size = lib_size[keep], stringsAsFactors = FALSE
  )
  rm(pb); gc(verbose = FALSE)
}

flag_all <- do.call(rbind, flag_rows)
write.table(flag_all, file.path(out_dir, "mesenchymal_flag_cells.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## --- Step 3: flagged fraction per cohort x condition -----------------------
flag_summary <- aggregate(flagged ~ cohort + condition, data = flag_all,
                           FUN = function(x) c(n_total = length(x), n_flagged = sum(x)))
flag_summary <- do.call(data.frame, flag_summary)
names(flag_summary) <- c("cohort", "condition", "n_total", "n_flagged")
flag_summary$pct_flagged <- 100 * flag_summary$n_flagged / flag_summary$n_total
flag_summary <- flag_summary[order(flag_summary$cohort, flag_summary$condition), ]
write.table(flag_summary, file.path(out_dir, "mesenchymal_flag_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== Flagged fraction per cohort x condition ===\n")
print(flag_summary, row.names = FALSE)

## --- Step 4: clean pseudobulk sample table + checkpoint --------------------
sample_table_clean <- do.call(rbind, sample_rows_clean)
sample_table_clean$sample_name <- paste(sample_table_clean$cohort, sample_table_clean$donor_id, sample_table_clean$condition, sep = "_")
for (co in names(pb_counts_clean)) {
  stopifnot(identical(colnames(pb_counts_clean[[co]]),
                       paste(sample_table_clean$donor_id[sample_table_clean$cohort == co],
                             sample_table_clean$condition[sample_table_clean$cohort == co], sep = "\r")))
  colnames(pb_counts_clean[[co]]) <- sample_table_clean$sample_name[sample_table_clean$cohort == co]
}
sample_table_clean <- sample_table_clean[order(sample_table_clean$cohort, sample_table_clean$donor_id, sample_table_clean$condition), ]
write.table(sample_table_clean, file.path(out_dir, "pseudobulk_clean_sample_table.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== Clean (unflagged) pseudobulk sample table ===\n")
print(sample_table_clean, row.names = FALSE)
cat(sprintf("\nTotal clean samples: %d (A=%d, B=%d, C=%d)\n", nrow(sample_table_clean),
            sum(sample_table_clean$cohort == "A"), sum(sample_table_clean$cohort == "B"), sum(sample_table_clean$cohort == "C")))

## Paired-donor recount, before (all-Endothelial, from the earlier checkpoint) vs after
orig_samples <- readRDS(file.path(SCRATCH_DIR, "pseudobulk_endothelial.rds"))$sample_table
count_paired <- function(st) {
  key <- paste(st$cohort, st$donor_id)
  tab <- table(key)
  by_cohort <- sapply(c("A", "B", "C"), function(co) sum(tab[grepl(paste0("^", co, " "), names(tab))] == 2))
  c(total = sum(tab == 2), by_cohort)
}
cat("\nPaired donors BEFORE (all-Endothelial):\n"); print(count_paired(orig_samples))
cat("Paired donors AFTER (clean, unflagged-Endothelial):\n"); print(count_paired(sample_table_clean))

saveRDS(list(counts = pb_counts_clean, sample_table = sample_table_clean),
        file.path(SCRATCH_DIR, "pseudobulk_endothelial_clean.rds"))
cat("\nSaved:", file.path(SCRATCH_DIR, "pseudobulk_endothelial_clean.rds"), "\n")
cat("Done.\n")
