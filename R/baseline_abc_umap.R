## Unintegrated A+B+C baseline: ScaleData + RunPCA + RunUMAP + DimPlot by
## cohort, mirroring baseline_scale_pca.R/baseline_umap.R but for all three
## cohorts. merged_ABC.rds is already subsetted to the 2000 shared
## integration features (done in prepare_integration_abc.R, before merging),
## so no further subsetting is needed here.
##
## Input:  merged_ABC.rds
## Output: results/figures/integration/ABC_umap_unintegrated_by_cohort.png
##         results/figures/integration/ABC_elbow_plot.png
##
## Run: Rscript R/baseline_abc_umap.R

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
library(ggplot2)

n_dims <- 10 # matches every other embedding in this project

cat("Loading merged_ABC.rds...\n")
merged_abc <- readRDS(file.path(SCRATCH_DIR, "merged_ABC.rds"))
cat(sprintf("  %d genes x %d cells\n", nrow(merged_abc), ncol(merged_abc)))

cat("ScaleData + RunPCA...\n")
merged_abc <- ScaleData(merged_abc, verbose = FALSE)
merged_abc <- RunPCA(merged_abc, npcs = 30, verbose = FALSE)

fig_dir <- file.path(repo, "results", "figures", "integration")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
p_elbow <- ElbowPlot(merged_abc, ndims = 30) + ggtitle("Unintegrated [A+B+C]: PCA elbow plot")
ggsave(file.path(fig_dir, "ABC_elbow_plot.png"), p_elbow, width = 7, height = 5, dpi = 150, bg = "white")

cat(sprintf("RunUMAP(dims = 1:%d)...\n", n_dims))
merged_abc <- RunUMAP(merged_abc, dims = 1:n_dims, verbose = FALSE)

abc_colors <- c(A = "#F8766D", B = "#00BFC4", C = "#7CAE00")
p <- DimPlot(merged_abc, reduction = "umap", group.by = "cohort") +
  scale_color_manual(values = abc_colors) +
  ggtitle("Unintegrated [A+B+C]: UMAP by cohort")
ggsave(file.path(fig_dir, "ABC_umap_unintegrated_by_cohort.png"), p, width = 8, height = 6, dpi = 150, bg = "white")

out_path <- file.path(SCRATCH_DIR, "merged_ABC_umap.rds")
saveRDS(merged_abc, out_path)
cat("\nSaved:", out_path, "\n")
cat("Wrote results/figures/integration/ABC_umap_unintegrated_by_cohort.png and ABC_elbow_plot.png\n")
