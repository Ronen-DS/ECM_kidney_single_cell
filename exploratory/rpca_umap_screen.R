## Classic Seurat anchor-based integration, step 3 (fast screening): on the
## "integrated" assay, ScaleData + RunPCA (a fresh PCA IS required here,
## unlike the Harmony path - the integrated assay's values themselves changed,
## so the old PCA is no longer valid) + RunUMAP, straight to a by-cohort
## DimPlot - skipping the neighbor-graph/resolution-sweep step for a fast
## answer on whether this approach mixes the cohorts better than Harmony did.
##
## Input:  integrated_AB_rpca.rds
## Output: merged_AB_rpca_umap.rds
##         results/figures/integration/AB_umap_rpca_by_cohort.png
##         results/figures/integration/AB_elbow_plot_rpca.png
##
## Run: Rscript R/rpca_umap_screen.R

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

cat("Loading integrated_AB_rpca.rds...\n")
integrated <- readRDS(file.path(SCRATCH_DIR, "integrated_AB_rpca.rds"))
DefaultAssay(integrated) <- "integrated"

cat("ScaleData + RunPCA on integrated assay...\n")
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)

fig_dir <- file.path(repo, "results", "figures", "integration")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
p_elbow <- ElbowPlot(integrated, ndims = 30) + ggtitle("Integrated [A+B] (RPCA anchors): elbow plot")
ggsave(file.path(fig_dir, "AB_elbow_plot_rpca.png"), p_elbow, width = 7, height = 5, dpi = 150, bg = "white")

n_dims <- 10 # reused from the unintegrated baseline / Harmony attempts, for consistency
cat(sprintf("RunUMAP(dims = 1:%d)...\n", n_dims))
integrated <- RunUMAP(integrated, dims = 1:n_dims, verbose = FALSE)

p <- DimPlot(integrated, reduction = "umap", group.by = "cohort") +
  ggtitle("Integrated [A+B] (RPCA anchors): UMAP by cohort")
ggsave(file.path(fig_dir, "AB_umap_rpca_by_cohort.png"), p, width = 7, height = 6, dpi = 150, bg = "white")

out_path <- file.path(SCRATCH_DIR, "merged_AB_rpca_umap.rds")
saveRDS(integrated, out_path)
cat("\nSaved:", out_path, "\n")
cat("Wrote results/figures/integration/AB_umap_rpca_by_cohort.png and AB_elbow_plot_rpca.png\n")
