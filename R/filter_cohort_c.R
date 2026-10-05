## Cohort C post-QC filtering: specimen-level minimum-cell filter, then a
## per-cell nFeature_RNA/non-MT-UMIs/percent_mt filter, applied to the
## seurat_C_raw.rds checkpoint. The raw nCount floor is replaced by a non-MT
## UMI floor (nCount * (1 - percent_mt/100) > 150) - the actual signal left
## once mitochondrial reads are set aside, which subsumes an nCount > 150
## floor automatically. The %MT cap comes from binning %MT and checking where
## median non-MT UMIs drops toward that same 150 floor (see
## derive_c_mt_threshold.R): every bin up to 90-95% still had a comfortable
## median (>=200), only 95-100% dropped close to it (160 overall, 153 for
## Adjacent) - so the cap is set at that bin edge, 95%, rather than the much
## more aggressive 20% used previously.
##   n_cells per specimen > 20
##   nFeature_RNA        > 100 & < 400
##   non-MT UMIs         > 150
##   percent_mt          < 95%
##
## Run: Rscript R/filter_cohort_c.R  (needs seurat_C_raw.rds - see R/build_seurat_c.R)

.args <- commandArgs(trailingOnly = FALSE)
.script_path <- sub("--file=", "", .args[grep("--file=", .args)])
repo <- normalizePath(file.path(dirname(.script_path), ".."))
readRenviron(file.path(repo, ".Renviron"))
source(file.path(repo, "config.R"))
source(file.path(repo, "R", "qc_utils.R"))
library(Seurat)
library(ggplot2)

cache_path <- file.path(SCRATCH_DIR, "seurat_C_raw.rds")
if (!file.exists(cache_path)) {
  stop("seurat_C_raw.rds not found - run R/build_seurat_c.R first.")
}
cat("Loading cached Seurat object...\n")
obj <- readRDS(cache_path)

result <- filter_seurat_cohort(
  obj,
  min_specimen_cells = 20,
  n_feature_min = 100, n_feature_max = 400,
  n_count_min = -Inf, n_count_max = Inf, # superseded by non_mt_umis_min below
  non_mt_umis_min = 150,
  percent_mt_max = 95
)

cat(sprintf(
  "Specimens: %d -> %d (dropped: %s)\n",
  result$n_specimens_before, result$n_specimens_after,
  if (length(result$dropped_specimens) > 0) paste(result$dropped_specimens, collapse = ", ") else "none"
))
cat(sprintf("Cells: %d -> %d (%.1f%% retained)\n",
            result$n_cells_before, result$n_cells_after,
            100 * result$n_cells_after / result$n_cells_before))

filtered <- result$obj
out_path <- file.path(SCRATCH_DIR, "seurat_C_filtered.rds")
saveRDS(filtered, out_path)
cat("Saved:", out_path, "\n")

## --- Post-filter QC table + figures, mirroring qc_cohort_c.R ---------------
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
write.table(per_specimen, file.path(tables_dir, "filter_cohort_c_per_specimen.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

fig_dir <- file.path(repo, "results", "figures", "post_qc", "C")
write_standard_qc_figures(qc, per_specimen, fig_dir,
                           cohort_label = "Cohort C (post-filter)",
                           file_prefix = "qc_cohort_c")

cat("\nWrote 6 post-filter figures and 1 table to results/figures/post_qc/C and results/tables/post_qc.\n")
