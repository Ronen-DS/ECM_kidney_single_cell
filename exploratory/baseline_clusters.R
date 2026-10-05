## Unintegrated baseline, part 2: FindNeighbors + a resolution sweep, using
## n_dims chosen from the elbow plot (results/figures/integration/AB_elbow_plot.png).
## Prints/saves cluster-count tables for each resolution so a resolution can
## be picked by inspection (not chosen automatically here) before running
## RunUMAP in the next step.
##
## Input:  merged_AB_pca.rds
## Output: merged_AB_neighbors.rds
##         results/tables/integration/AB_resolution_sweep.tsv
##
## Run: Rscript R/baseline_clusters.R

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

n_dims_unint <- 10 # chosen from AB_elbow_plot.png

cat("Loading merged_AB_pca.rds...\n")
merged_ab <- readRDS(file.path(SCRATCH_DIR, "merged_AB_pca.rds"))

cat(sprintf("FindNeighbors(dims = 1:%d)...\n", n_dims_unint))
merged_ab <- FindNeighbors(merged_ab, dims = 1:n_dims_unint, verbose = FALSE)

res_vals <- c(0.3, 0.5, 0.7, 1.0)
sweep_rows <- list()
for (r in res_vals) {
  merged_ab <- FindClusters(merged_ab, resolution = r, verbose = FALSE)
  counts <- table(Idents(merged_ab))
  cat(sprintf("\nUnintegrated - Resolution: %s (%d clusters)\n", r, length(counts)))
  print(counts)
  sweep_rows[[length(sweep_rows) + 1]] <- data.frame(
    resolution = r, n_clusters = length(counts),
    min_cluster_size = min(counts), max_cluster_size = max(counts)
  )
}
sweep_tbl <- do.call(rbind, sweep_rows)

tables_dir <- file.path(repo, "results", "tables", "integration")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
write.table(sweep_tbl, file.path(tables_dir, "AB_resolution_sweep.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

out_path <- file.path(SCRATCH_DIR, "merged_AB_neighbors.rds")
saveRDS(merged_ab, out_path)
cat("\nSaved:", out_path, "\n")
cat("\nResolution sweep summary:\n")
print(sweep_tbl, row.names = FALSE)
cat("\nPick a resolution for RunUMAP/final clustering based on this table.\n")
