## RPCA-integrated clustering + final plots, mirroring baseline_clusters.R /
## baseline_umap.R and harmony_clusters.R / harmony_umap.R but on the RPCA-
## integrated object. UMAP was already computed in rpca_umap_screen.R (on the
## integrated assay's own "pca" reduction, dims=1:10) so it's reused directly
## here rather than recomputed - only FindNeighbors/FindClusters are new.
##
## Resolution sweep limited to 0.5 and 1.0 only (per user request, to save
## time) rather than the full 0.3/0.5/0.7/1.0 sweep used for the unintegrated
## baseline and Harmony.
##
## Input:  merged_AB_rpca_umap.rds
## Output: merged_AB_rpca_final.rds
##         results/tables/integration/AB_resolution_sweep_rpca.tsv
##         results/figures/integration/AB_umap_rpca_res0.5.png
##         results/figures/integration/AB_umap_rpca_res1.png
##
## Run: Rscript R/rpca_clusters_umap.R

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

n_dims <- 10 # matches the unintegrated baseline / Harmony attempts

cat("Loading merged_AB_rpca_umap.rds...\n")
merged_ab <- readRDS(file.path(SCRATCH_DIR, "merged_AB_rpca_umap.rds"))

cat(sprintf("FindNeighbors(reduction = 'pca', dims = 1:%d)...\n", n_dims))
merged_ab <- FindNeighbors(merged_ab, reduction = "pca", dims = 1:n_dims,
                            graph.name = c("rpca_nn", "rpca_snn"), verbose = FALSE)

res_vals <- c(0.5, 1.0)
sweep_rows <- list()
for (r in res_vals) {
  merged_ab <- FindClusters(merged_ab, resolution = r, graph.name = "rpca_snn", verbose = FALSE)
  col <- paste0("rpca_snn_res.", r)
  counts <- table(merged_ab@meta.data[[col]])
  cat(sprintf("\nIntegrated (RPCA) - Resolution: %s (%d clusters)\n", r, length(counts)))
  print(counts)
  sweep_rows[[length(sweep_rows) + 1]] <- data.frame(
    resolution = r, n_clusters = length(counts),
    min_cluster_size = min(counts), max_cluster_size = max(counts)
  )
}
sweep_tbl <- do.call(rbind, sweep_rows)

tables_dir <- file.path(repo, "results", "tables", "integration")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
write.table(sweep_tbl, file.path(tables_dir, "AB_resolution_sweep_rpca.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

## --- Plots (reusing the UMAP already computed in rpca_umap_screen.R) -------
fig_dir <- file.path(repo, "results", "figures", "integration")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

Idents(merged_ab) <- merged_ab$rpca_snn_res.0.5
p_res05 <- DimPlot(merged_ab, reduction = "umap", label = TRUE) +
  ggtitle("Integrated [A+B] (RPCA anchors): UMAP, clusters at resolution 0.5")
ggsave(file.path(fig_dir, "AB_umap_rpca_res0.5.png"), p_res05, width = 7, height = 6, dpi = 150, bg = "white")

Idents(merged_ab) <- merged_ab$rpca_snn_res.1
p_res1 <- DimPlot(merged_ab, reduction = "umap", label = TRUE) +
  ggtitle("Integrated [A+B] (RPCA anchors): UMAP, clusters at resolution 1.0")
ggsave(file.path(fig_dir, "AB_umap_rpca_res1.png"), p_res1, width = 7, height = 6, dpi = 150, bg = "white")

Idents(merged_ab) <- merged_ab$rpca_snn_res.0.5
merged_ab$seurat_clusters_rpca <- merged_ab$rpca_snn_res.0.5

out_path <- file.path(SCRATCH_DIR, "merged_AB_rpca_final.rds")
saveRDS(merged_ab, out_path)
cat("\nSaved:", out_path, "\n")
cat("\nResolution sweep summary (RPCA):\n")
print(sweep_tbl, row.names = FALSE)
cat("\nWrote 2 UMAP figures to results/figures/integration/\n")
