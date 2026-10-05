## Derives cohort C's %MT cutoff empirically, tied directly to the nCount
## floor already chosen, instead of an arbitrary round-number cap:
##   1. Apply the nFeature (100-400) and nCount (>150) filter only.
##   2. On the cells that remain, compute each cell's non-MT UMI count -
##      nCount * (1 - percent_mt/100) - the signal actually left for calling
##      identity once mitochondrial reads are set aside.
##   3. Bin by percent_mt (0-20, 20-50, 50-70, 70-80, 80-90, 90-95, 95-100)
##      and report median non-MT UMIs and n_cells per bin, split by
##      condition (Tumor/Adjacent) as a robustness check that the pattern
##      isn't driven by one condition alone.
##   4. The %MT cap is the lower edge of the first bin (in increasing %MT
##      order) whose OVERALL median non-MT UMIs drops below 150 (the same
##      floor already used for nCount) - i.e. cells above this %MT typically
##      have fewer nuclear UMIs than the minimum required of any cell.
##
## This is a diagnostic/derivation step only - it does not write a new
## filtered checkpoint. See filter_cohort_c.R for where the derived cutoff
## gets applied.
##
## Input:  seurat_C_raw.rds
## Output: results/tables/post_qc/cohort_c_mt_binning.tsv
##
## Run: Rscript R/derive_c_mt_threshold.R

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

N_FEATURE_MIN <- 100
N_FEATURE_MAX <- 400
N_COUNT_MIN <- 150

cat("Loading seurat_C_raw.rds...\n")
obj <- readRDS(file.path(SCRATCH_DIR, "seurat_C_raw.rds"))
meta <- obj@meta.data
cat(sprintf("Raw cells: %d\n", nrow(meta)))

pass <- meta$nFeature_RNA > N_FEATURE_MIN & meta$nFeature_RNA < N_FEATURE_MAX & meta$nCount_RNA > N_COUNT_MIN
qc <- meta[pass, ]
cat(sprintf("After nFeature (%d-%d) + nCount (>%d) filter: %d cells (%.1f%%)\n",
            N_FEATURE_MIN, N_FEATURE_MAX, N_COUNT_MIN, nrow(qc), 100 * nrow(qc) / nrow(meta)))

qc$non_mt_umis <- qc$nCount_RNA * (1 - qc$percent_mt / 100)

mt_breaks <- c(0, 20, 50, 70, 80, 90, 95, 100)
mt_labels <- c("0-20", "20-50", "50-70", "70-80", "80-90", "90-95", "95-100")
qc$mt_bin <- cut(qc$percent_mt, breaks = mt_breaks, labels = mt_labels, include.lowest = TRUE, right = TRUE)

## Overall (drives the cutoff rule) ------------------------------------------
overall <- aggregate(non_mt_umis ~ mt_bin, data = qc, FUN = median)
overall_n <- as.data.frame(table(mt_bin = qc$mt_bin))
overall <- merge(overall, overall_n, by = "mt_bin")
colnames(overall) <- c("mt_bin", "median_non_mt_umis", "n_cells")
overall <- overall[match(mt_labels, overall$mt_bin), ] # keep bin order

cat("\n=== Overall (drives the cutoff) ===\n")
print(overall, row.names = FALSE)

## Per-condition (reported for robustness, per the assignment) ---------------
by_condition <- aggregate(non_mt_umis ~ mt_bin + condition, data = qc, FUN = median)
n_by_condition <- as.data.frame(table(mt_bin = qc$mt_bin, condition = qc$condition))
colnames(n_by_condition) <- c("mt_bin", "condition", "n_cells")
by_condition <- merge(by_condition, n_by_condition, by = c("mt_bin", "condition"))
colnames(by_condition)[colnames(by_condition) == "non_mt_umis"] <- "median_non_mt_umis"
by_condition <- by_condition[order(match(by_condition$mt_bin, mt_labels), by_condition$condition), ]

cat("\n=== Split by condition ===\n")
print(by_condition, row.names = FALSE)

## Derive the cutoff -----------------------------------------------------------
below_floor <- which(overall$median_non_mt_umis < N_COUNT_MIN)
if (length(below_floor) == 0) {
  cat(sprintf("\nNo bin's median non-MT UMIs drops below %d - no data-driven cap found in this range.\n", N_COUNT_MIN))
} else {
  first_bin_idx <- below_floor[1]
  mt_cap <- mt_breaks[first_bin_idx] # lower edge of that bin
  cat(sprintf(
    "\nFirst bin below the nCount floor (%d): '%s' (median non-MT UMIs = %.1f, n = %d)\n",
    N_COUNT_MIN, overall$mt_bin[first_bin_idx], overall$median_non_mt_umis[first_bin_idx], overall$n_cells[first_bin_idx]
  ))
  cat(sprintf("Derived %%MT cap: percent_mt < %d\n", mt_cap))
  cat("Justification: cells above this %MT typically have fewer nuclear UMIs than the minimum we require of any cell.\n")
}

tables_dir <- file.path(repo, "results", "tables", "post_qc")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
write.table(overall, file.path(tables_dir, "cohort_c_mt_binning.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(by_condition, file.path(tables_dir, "cohort_c_mt_binning_by_condition.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nWrote results/tables/post_qc/cohort_c_mt_binning.tsv and cohort_c_mt_binning_by_condition.tsv\n")
