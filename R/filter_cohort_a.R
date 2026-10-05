## Cohort A post-QC filtering: specimen-level minimum-cell filter, then a
## per-nucleus nFeature_RNA/nCount_RNA/percent_mt filter, applied to the
## seurat_A_raw.rds checkpoint. Thresholds are user-specified (not re-derived
## here):
##   n_cells per specimen > 20
##   nFeature_RNA > 500 & < 4500
##   nCount_RNA   > 1000
##   percent_mt   < 2%
##
## Run: Rscript R/filter_cohort_a.R  (needs seurat_A_raw.rds - see R/build_seurat_a.R)

.args <- commandArgs(trailingOnly = FALSE)
.script_path <- sub("--file=", "", .args[grep("--file=", .args)])
repo <- normalizePath(file.path(dirname(.script_path), ".."))
readRenviron(file.path(repo, ".Renviron"))
source(file.path(repo, "config.R"))
source(file.path(repo, "R", "qc_utils.R"))
library(Seurat)
library(ggplot2)

cache_path <- file.path(SCRATCH_DIR, "seurat_A_raw.rds")
if (!file.exists(cache_path)) {
  stop("seurat_A_raw.rds not found - run R/build_seurat_a.R first.")
}
cat("Loading cached Seurat object...\n")
obj <- readRDS(cache_path)

result <- filter_seurat_cohort(
  obj,
  min_specimen_cells = 20,
  n_feature_min = 500, n_feature_max = 4500,
  n_count_min = 1000, n_count_max = Inf,
  percent_mt_max = 2
)

cat(sprintf(
  "Specimens: %d -> %d (dropped: %s)\n",
  result$n_specimens_before, result$n_specimens_after,
  if (length(result$dropped_specimens) > 0) paste(result$dropped_specimens, collapse = ", ") else "none"
))
cat(sprintf("Nuclei: %d -> %d (%.1f%% retained)\n",
            result$n_cells_before, result$n_cells_after,
            100 * result$n_cells_after / result$n_cells_before))

filtered <- result$obj
out_path <- file.path(SCRATCH_DIR, "seurat_A_filtered.rds")
saveRDS(filtered, out_path)
cat("Saved:", out_path, "\n")

## --- Post-filter QC table + figures, mirroring qc_cohort_a.R ---------------
qc <- filtered@meta.data
qc$n_count <- qc$nCount_RNA
qc$n_feature <- qc$nFeature_RNA
qc$sample_id <- as.character(qc$sample_id)

per_specimen <- aggregate(
  cbind(n_count, n_feature, percent_mt) ~ sample_id + donor_id + condition,
  data = qc, FUN = median
)
n_cells_tab <- as.data.frame(table(sample_id = qc$sample_id))
colnames(n_cells_tab) <- c("sample_id", "n_cells")
per_specimen <- merge(per_specimen, n_cells_tab, by = "sample_id")
colnames(per_specimen)[colnames(per_specimen) %in% c("n_count", "n_feature", "percent_mt")] <-
  c("median_n_count", "median_n_feature", "median_percent_mt")
per_specimen <- per_specimen[order(per_specimen$n_cells), ]

tables_dir <- file.path(repo, "results", "tables", "post_qc")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
write.table(per_specimen, file.path(tables_dir, "filter_cohort_a_per_specimen.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

fig_dir <- file.path(repo, "results", "figures", "post_qc", "A")
write_standard_qc_figures(qc, per_specimen, fig_dir,
                           cohort_label = "Cohort A (post-filter)",
                           file_prefix = "qc_cohort_a")

cat("\nWrote 6 post-filter figures and 1 table to results/figures/post_qc/A and results/tables/post_qc.\n")
