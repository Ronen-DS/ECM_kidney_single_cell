## Unintegrated baseline, part 1: ScaleData + RunPCA + ElbowPlot on the
## merged [A+B] object, using the SelectIntegrationFeatures() gene set
## already computed in prepare_integration.R (not each cohort's own
## independently-chosen HVGs - the whole point of that step was to pick a
## feature set valid across both cohorts for exactly this embedding).
##
## Layers are left split by cohort (Seurat v5) rather than joined - this
## matches the standard Seurat v5 "before integration" workflow, and the
## same "pca" reduction computed here is reused as the input Harmony
## integrates FROM in the next step (IntegrateLayers() corrects an existing
## PCA rather than requiring a second one on a separate "integrated" assay).
##
## Subsets to the 2000 integration features BEFORE scaling (not the full
## 45,799-gene union) - this cuts the in-memory matrix size by ~20x and is
## what caused the earlier near-OOM to actually be avoidable: PCA only ever
## uses this feature subset anyway, so there is no reason to carry the full
## gene set through this step.
##
## Input:  merged_AB.rds, results/tables/integration/integration_features_AB.txt
## Output: merged_AB_pca.rds
##         results/figures/integration/AB_elbow_plot.png
##
## Run: Rscript R/baseline_scale_pca.R

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

cat("Loading merged_AB.rds...\n")
merged_ab <- readRDS(file.path(SCRATCH_DIR, "merged_AB.rds"))
cat(sprintf("  Loaded: %d genes x %d cells\n", nrow(merged_ab), ncol(merged_ab)))

features_ab <- readLines(file.path(repo, "results", "tables", "integration", "integration_features_AB.txt"))
cat(sprintf("Subsetting to %d integration features...\n", length(features_ab)))
merged_ab <- subset(merged_ab, features = features_ab)
cat(sprintf("  Subsetted: %d genes x %d cells\n", nrow(merged_ab), ncol(merged_ab)))

cat("ScaleData...\n")
merged_ab <- ScaleData(merged_ab, features = features_ab, verbose = FALSE)

cat("RunPCA...\n")
merged_ab <- RunPCA(merged_ab, features = features_ab, npcs = 50, verbose = FALSE)

fig_dir <- file.path(repo, "results", "figures", "integration")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
p <- ElbowPlot(merged_ab, ndims = 50) + ggtitle("Unintegrated [A+B]: PCA elbow plot")
ggsave(file.path(fig_dir, "AB_elbow_plot.png"), p, width = 7, height = 5, dpi = 150, bg = "white")

out_path <- file.path(SCRATCH_DIR, "merged_AB_pca.rds")
saveRDS(merged_ab, out_path)
cat("\nSaved:", out_path, "\n")
cat("Wrote results/figures/integration/AB_elbow_plot.png - inspect to choose n_dims for FindNeighbors/RunUMAP.\n")
