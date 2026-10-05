## Adds resolution 0.3 to the A+B+C RPCA resolution sweep. FindNeighbors'
## graph (rpca_snn) is already stored in merged_ABC_rpca_final.rds, so only
## FindClusters (cheap) needs to run - no need to redo the expensive
## neighbor search from rpca_clusters_umap_abc.R.
##
## Input:  merged_ABC_rpca_final.rds
## Output: merged_ABC_rpca_final.rds (updated in place, adds resolution 0.3)
##         results/tables/integration/ABC_resolution_sweep_rpca.tsv (updated)
##         results/figures/integration/ABC_umap_rpca_res0.3.png
##
## Run: Rscript R/rpca_cluster_res03_abc.R

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

cat("Loading merged_ABC_rpca_final.rds...\n")
merged_abc <- readRDS(file.path(SCRATCH_DIR, "merged_ABC_rpca_final.rds"))

cat("FindClusters(resolution = 0.3, graph.name = 'rpca_snn')...\n")
merged_abc <- FindClusters(merged_abc, resolution = 0.3, graph.name = "rpca_snn", verbose = FALSE)
counts <- table(merged_abc$rpca_snn_res.0.3)
cat(sprintf("\nIntegrated (RPCA, A+B+C) - Resolution: 0.3 (%d clusters)\n", length(counts)))
print(counts)

fig_dir <- file.path(repo, "results", "figures", "integration")
Idents(merged_abc) <- merged_abc$rpca_snn_res.0.3
p <- DimPlot(merged_abc, reduction = "umap", label = TRUE) +
  ggtitle("Integrated [A+B+C] (RPCA anchors): UMAP, clusters at resolution 0.3")
ggsave(file.path(fig_dir, "ABC_umap_rpca_res0.3.png"), p, width = 7, height = 6, dpi = 150, bg = "white")

## --- Update the resolution sweep table (adds 0.3, keeps existing 0.5/1.0 rows) ---
tables_dir <- file.path(repo, "results", "tables", "integration")
sweep_path <- file.path(tables_dir, "ABC_resolution_sweep_rpca.tsv")
existing <- read.delim(sweep_path, stringsAsFactors = FALSE)
new_row <- data.frame(resolution = 0.3, n_clusters = length(counts),
                       min_cluster_size = min(counts), max_cluster_size = max(counts))
sweep_tbl <- rbind(new_row, existing)
sweep_tbl <- sweep_tbl[order(sweep_tbl$resolution), ]
write.table(sweep_tbl, sweep_path, sep = "\t", row.names = FALSE, quote = FALSE)

Idents(merged_abc) <- merged_abc$rpca_snn_res.0.5 # restore the previously-designated default
out_path <- file.path(SCRATCH_DIR, "merged_ABC_rpca_final.rds")
saveRDS(merged_abc, out_path)
cat("\nSaved:", out_path, "\n")
cat("\nUpdated resolution sweep summary (RPCA, A+B+C):\n")
print(sweep_tbl, row.names = FALSE)
cat("\nWrote results/figures/integration/ABC_umap_rpca_res0.3.png\n")
