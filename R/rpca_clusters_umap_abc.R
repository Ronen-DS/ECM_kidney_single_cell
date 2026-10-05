## RPCA-integrated A+B+C clustering + final plots, mirroring
## rpca_clusters_umap.R (the A+B version) but on the joint A+B+C integration.
## UMAP was already computed in rpca_umap_only_abc.R (on the integrated
## assay's own "pca" reduction, dims=1:10) so it's reused directly - only
## FindNeighbors/FindClusters are new. Resolutions limited to 0.5 and 1.0,
## matching the A+B version.
##
## Input:  merged_ABC_rpca_umap.rds
## Output: merged_ABC_rpca_final.rds
##         results/tables/integration/ABC_resolution_sweep_rpca.tsv
##         results/figures/integration/ABC_umap_rpca_res0.5.png
##         results/figures/integration/ABC_umap_rpca_res1.png
##
## Run: Rscript R/rpca_clusters_umap_abc.R

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

n_dims <- 10

cat("Loading merged_ABC_rpca_umap.rds...\n")
merged_abc <- readRDS(file.path(SCRATCH_DIR, "merged_ABC_rpca_umap.rds"))

cat(sprintf("FindNeighbors(reduction = 'pca', dims = 1:%d)...\n", n_dims))
merged_abc <- FindNeighbors(merged_abc, reduction = "pca", dims = 1:n_dims,
                             graph.name = c("rpca_nn", "rpca_snn"), verbose = FALSE)

res_vals <- c(0.5, 1.0)
sweep_rows <- list()
for (r in res_vals) {
  merged_abc <- FindClusters(merged_abc, resolution = r, graph.name = "rpca_snn", verbose = FALSE)
  col <- paste0("rpca_snn_res.", r)
  counts <- table(merged_abc@meta.data[[col]])
  cat(sprintf("\nIntegrated (RPCA, A+B+C) - Resolution: %s (%d clusters)\n", r, length(counts)))
  print(counts)
  sweep_rows[[length(sweep_rows) + 1]] <- data.frame(
    resolution = r, n_clusters = length(counts),
    min_cluster_size = min(counts), max_cluster_size = max(counts)
  )
}
sweep_tbl <- do.call(rbind, sweep_rows)

tables_dir <- file.path(repo, "results", "tables", "integration")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
write.table(sweep_tbl, file.path(tables_dir, "ABC_resolution_sweep_rpca.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

## --- Plots (reusing the UMAP already computed) -------------------------------
fig_dir <- file.path(repo, "results", "figures", "integration")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

Idents(merged_abc) <- merged_abc$rpca_snn_res.0.5
p_res05 <- DimPlot(merged_abc, reduction = "umap", label = TRUE) +
  ggtitle("Integrated [A+B+C] (RPCA anchors): UMAP, clusters at resolution 0.5")
ggsave(file.path(fig_dir, "ABC_umap_rpca_res0.5.png"), p_res05, width = 7, height = 6, dpi = 150, bg = "white")

Idents(merged_abc) <- merged_abc$rpca_snn_res.1
p_res1 <- DimPlot(merged_abc, reduction = "umap", label = TRUE) +
  ggtitle("Integrated [A+B+C] (RPCA anchors): UMAP, clusters at resolution 1.0")
ggsave(file.path(fig_dir, "ABC_umap_rpca_res1.png"), p_res1, width = 7, height = 6, dpi = 150, bg = "white")

Idents(merged_abc) <- merged_abc$rpca_snn_res.0.5
merged_abc$seurat_clusters_rpca <- merged_abc$rpca_snn_res.0.5

out_path <- file.path(SCRATCH_DIR, "merged_ABC_rpca_final.rds")
saveRDS(merged_abc, out_path)
cat("\nSaved:", out_path, "\n")
cat("\nResolution sweep summary (RPCA, A+B+C):\n")
print(sweep_tbl, row.names = FALSE)
cat("\nWrote 2 UMAP figures to results/figures/integration/\n")
