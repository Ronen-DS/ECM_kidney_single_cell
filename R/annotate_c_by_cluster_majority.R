## Final cell annotation for the A+B+C integrated UMAP:
##   - A and B cells: harmonized major_population directly from
##     label_mapping.tsv (via add_major_population()).
##   - C cells: assigned the majority A/B population within their
##     resolution-1 cluster of the RPCA-integrated graph (finer clusters are
##     more homogeneous, maximizing assignment purity). Purity is computed
##     on A/B cells only (C cells never vote on their own label, and never
##     contribute to another C cell's label).
##     A cluster assigns a label to its C cells only if:
##       - the cluster has >= 50 A/B cells, AND
##       - the majority A/B population makes up >= 80% of those A/B cells, AND
##       - that majority population is not itself "Unresolved".
##     Otherwise (too few A/B cells, insufficient purity, or an Unresolved
##     majority), that cluster's C cells are marked "Unresolved".
##
## Input:  merged_ABC_rpca_final.rds, metadata/label_mapping.tsv
## Output: merged_ABC_annotated.rds
##         results/tables/integration/ABC_cluster_majority_assignment.tsv
##         results/figures/integration/ABC_umap_final_annotation.png
##
## Run: Rscript R/annotate_c_by_cluster_majority.R

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
library(ggplot2)

MIN_AB_CELLS <- 50
MIN_PURITY_PCT <- 80

label_mapping <- read.delim(LABEL_MAPPING_PATH, stringsAsFactors = FALSE)

cat("Loading merged_ABC_rpca_final.rds...\n")
obj <- readRDS(file.path(SCRATCH_DIR, "merged_ABC_rpca_final.rds"))
obj <- add_major_population(obj, label_mapping) # A/B get their harmonized population; C gets "Unresolved" here (unused for C below)

cluster <- obj$rpca_snn_res.1
is_ab <- obj$cohort %in% c("A", "B")
is_c <- obj$cohort == "C"

cluster_levels <- sort(unique(cluster))
assignment_rows <- lapply(cluster_levels, function(cl) {
  ab_pop <- obj$major_population[is_ab & cluster == cl]
  n_ab <- length(ab_pop)
  n_c <- sum(is_c & cluster == cl)

  if (n_ab == 0) {
    majority_pop <- NA_character_
    majority_pct <- NA_real_
  } else {
    tab <- table(ab_pop)
    majority_pop <- names(tab)[which.max(tab)]
    majority_pct <- 100 * max(tab) / n_ab
  }

  assigned <- !is.na(majority_pop) && n_ab >= MIN_AB_CELLS &&
    majority_pct >= MIN_PURITY_PCT && majority_pop != "Unresolved"
  c_label <- if (assigned) majority_pop else "Unresolved"

  data.frame(
    cluster = cl, n_ab = n_ab, n_c = n_c,
    majority_population = ifelse(is.na(majority_pop), "(no A/B cells)", majority_pop),
    majority_pct = round(majority_pct, 1),
    assigned = assigned,
    c_label = c_label
  )
})
assignment_tbl <- do.call(rbind, assignment_rows)
assignment_tbl <- assignment_tbl[order(-assignment_tbl$n_c), ]

tables_dir <- file.path(repo, "results", "tables", "integration")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
write.table(assignment_tbl, file.path(tables_dir, "ABC_cluster_majority_assignment.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nPer-cluster majority assignment (resolution 1):\n")
print(assignment_tbl, row.names = FALSE)
cat(sprintf("\nClusters assigned a real label: %d / %d\n", sum(assignment_tbl$assigned), nrow(assignment_tbl)))

## --- Build the final combined annotation -------------------------------------
c_label_by_cluster <- setNames(assignment_tbl$c_label, assignment_tbl$cluster)
final_annotation <- obj$major_population
final_annotation[is_c] <- c_label_by_cluster[as.character(cluster[is_c])]
obj$cell_annotation <- final_annotation

cat("\nFinal annotation counts (all cells):\n")
print(table(obj$cell_annotation))
cat("\nFinal annotation counts, cohort C only:\n")
print(table(obj$cell_annotation[is_c]))

## --- Figure --------------------------------------------------------------
populations <- sort(unique(label_mapping$population))
pop_colors <- setNames(scales::hue_pal()(length(populations)), populations)
pop_colors["Unresolved"] <- "grey80"

fig_dir <- file.path(repo, "results", "figures", "integration")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
p <- DimPlot(obj, reduction = "umap", group.by = "cell_annotation") +
  scale_color_manual(values = pop_colors) +
  ggtitle("Integrated [A+B+C]: final cell annotation (A/B harmonized; C via resolution-1 cluster majority)")
ggsave(file.path(fig_dir, "ABC_umap_final_annotation.png"), p, width = 8.5, height = 6, dpi = 150, bg = "white")

out_path <- file.path(SCRATCH_DIR, "merged_ABC_annotated.rds")
saveRDS(obj, out_path)
cat("\nSaved:", out_path, "\n")
cat("Wrote results/figures/integration/ABC_umap_final_annotation.png\n")
