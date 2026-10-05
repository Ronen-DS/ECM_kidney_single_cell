## Unintegrated baseline, part 3: RunUMAP + DimPlots. Final chosen resolution
## is 0.5, but the resolution sweep in baseline_clusters.R already computed
## and stored cluster assignments for ALL of 0.3/0.5/0.7/1.0 as separate
## meta.data columns (RNA_snn_res.<x>) - UMAP itself only depends on the PCA
## dims (not resolution), so it's run once and both a res=0.5 and a res=1.0
## coloring are plotted from the assignments already on disk, no re-clustering
## needed.
##
## Input:  merged_AB_neighbors.rds
## Output: merged_AB_umap.rds
##         results/figures/integration/AB_umap_unintegrated_by_cohort.png
##         results/figures/integration/AB_umap_unintegrated_res0.5.png
##         results/figures/integration/AB_umap_unintegrated_res1.png
##
## Run: Rscript R/baseline_umap.R

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

n_dims_unint <- 10 # matches baseline_clusters.R

cat("Loading merged_AB_neighbors.rds...\n")
merged_ab <- readRDS(file.path(SCRATCH_DIR, "merged_AB_neighbors.rds"))

cat(sprintf("RunUMAP(dims = 1:%d)...\n", n_dims_unint))
merged_ab <- RunUMAP(merged_ab, dims = 1:n_dims_unint, verbose = FALSE)

fig_dir <- file.path(repo, "results", "figures", "integration")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

p_cohort <- DimPlot(merged_ab, reduction = "umap", group.by = "cohort") +
  ggtitle("Unintegrated [A+B]: UMAP by cohort")
ggsave(file.path(fig_dir, "AB_umap_unintegrated_by_cohort.png"), p_cohort,
       width = 7, height = 6, dpi = 150, bg = "white")

Idents(merged_ab) <- merged_ab$RNA_snn_res.0.5
p_res05 <- DimPlot(merged_ab, reduction = "umap", label = TRUE) +
  ggtitle("Unintegrated [A+B]: UMAP, clusters at resolution 0.5")
ggsave(file.path(fig_dir, "AB_umap_unintegrated_res0.5.png"), p_res05,
       width = 7, height = 6, dpi = 150, bg = "white")

Idents(merged_ab) <- merged_ab$RNA_snn_res.1
p_res1 <- DimPlot(merged_ab, reduction = "umap", label = TRUE) +
  ggtitle("Unintegrated [A+B]: UMAP, clusters at resolution 1.0")
ggsave(file.path(fig_dir, "AB_umap_unintegrated_res1.png"), p_res1,
       width = 7, height = 6, dpi = 150, bg = "white")

## Chosen resolution (0.5) as the "final" Idents/seurat_clusters going forward
Idents(merged_ab) <- merged_ab$RNA_snn_res.0.5
merged_ab$seurat_clusters <- merged_ab$RNA_snn_res.0.5

out_path <- file.path(SCRATCH_DIR, "merged_AB_umap.rds")
saveRDS(merged_ab, out_path)
cat("\nSaved:", out_path, "\n")
cat("Wrote 3 UMAP figures to results/figures/integration/\n")
