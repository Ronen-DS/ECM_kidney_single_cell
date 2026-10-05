## Cross-tabulates resolution-0.3 clusters (A+B+C RPCA integration) against
## the harmonized cell population (A/B) / cohort C's own category, to check
## quantitatively whether 0.3's coarser clusters are actually dominated by
## one population each (clean correspondence) or still mix/split
## populations - not just an eyeballed impression from the UMAP colors.
##
## Input:  merged_ABC_rpca_final.rds, metadata/label_mapping.tsv
## Output: results/tables/integration/ABC_res0.3_population_crosstab.tsv
##         results/tables/integration/ABC_res0.3_cluster_purity.tsv
##
## Run: Rscript R/crosstab_res03_population_abc.R

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
source(file.path(repo, "R", "integration_utils.R"))
library(Seurat)

label_mapping <- read.delim(LABEL_MAPPING_PATH, stringsAsFactors = FALSE)

cat("Loading merged_ABC_rpca_final.rds...\n")
obj <- readRDS(file.path(SCRATCH_DIR, "merged_ABC_rpca_final.rds"))
obj <- add_major_population(obj, label_mapping)

is_c <- obj$cohort == "C"
plot_group <- obj$major_population
plot_group[is_c] <- "Cohort C"
obj$plot_group <- plot_group

cluster <- obj$rpca_snn_res.0.3
crosstab <- table(cluster = cluster, population = obj$plot_group)

tables_dir <- file.path(repo, "results", "tables", "integration")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
write.table(as.data.frame.matrix(crosstab), file.path(tables_dir, "ABC_res0.3_population_crosstab.tsv"),
            sep = "\t", quote = FALSE)

## Per-cluster purity: dominant category + its share of the cluster
purity_rows <- lapply(rownames(crosstab), function(cl) {
  row <- crosstab[cl, ]
  total <- sum(row)
  dominant <- names(row)[which.max(row)]
  data.frame(
    cluster = cl,
    n_cells = total,
    dominant_category = dominant,
    dominant_pct = round(100 * max(row) / total, 1)
  )
})
purity_tbl <- do.call(rbind, purity_rows)
purity_tbl <- purity_tbl[order(-purity_tbl$n_cells), ]
write.table(purity_tbl, file.path(tables_dir, "ABC_res0.3_cluster_purity.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nPer-cluster dominant category and purity (resolution 0.3):\n")
print(purity_tbl, row.names = FALSE)
cat(sprintf("\nMedian cluster purity: %.1f%%\n", median(purity_tbl$dominant_pct)))
cat(sprintf("Clusters with >=90%% purity: %d / %d\n", sum(purity_tbl$dominant_pct >= 90), nrow(purity_tbl)))

cat("\nFull crosstab:\n")
print(crosstab)
