## ECM panel expression across major cell populations, split into 3
## side-by-side panels (A, B, C) instead of pooled - the first step of the
## robustness check (does the pooled pattern from ecm_panel_by_population.R
## hold up within each cohort separately, or is it driven by just one).
##
## Relative expression (z) is computed WITHIN each cohort's own gene-wise
## distribution across its own populations - each panel is self-normalized,
## showing that cohort's own pattern, not a globally-pooled one. A
## cohort/population combination with zero cells (e.g. cohort C has no
## Stromal_fibro_myo cells at all - see the cluster-majority annotation
## step) is simply absent from that panel rather than shown as a fabricated
## zero.
##
## Input:  merged_ABC_annotated.rds, seurat_{A,B,C}_normalized.rds
## Output: results/tables/ecm_panel_by_population_per_cohort.tsv
##         results/figures/ecm_panel_by_population_per_cohort.png
##
## Run: Rscript R/ecm_panel_by_population_per_cohort.R

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
library(Matrix)
library(ggplot2)

ecm_genes <- c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1",
               "COL6A2", "FN1", "SPARC", "DCN", "LUM")

cat("Loading merged_ABC_annotated.rds for metadata only...\n")
obj_meta <- readRDS(file.path(SCRATCH_DIR, "merged_ABC_annotated.rds"))
meta <- data.frame(
  cell_key = colnames(obj_meta),
  cohort = obj_meta$cohort,
  cell_annotation = obj_meta$cell_annotation,
  stringsAsFactors = FALSE
)
rm(obj_meta)
gc(verbose = FALSE)
cat(sprintf("  %d cells\n", nrow(meta)))

ecm_frames <- list()
for (cohort in c("A", "B", "C")) {
  cat(sprintf("Cohort %s...\n", cohort))
  obj <- readRDS(file.path(SCRATCH_DIR, sprintf("seurat_%s_normalized.rds", cohort)))
  mat <- GetAssayData(obj, assay = "RNA", layer = "data")[ecm_genes, , drop = FALSE]
  df <- as.data.frame(as.matrix(t(mat)))
  df$cell_key <- colnames(obj)
  ecm_frames[[cohort]] <- df[, c("cell_key", ecm_genes)]
  rm(obj, mat)
  gc(verbose = FALSE)
}
ecm_df <- do.call(rbind, ecm_frames)
combined <- merge(meta, ecm_df, by = "cell_key")
cat(sprintf("Combined: %d cells\n", nrow(combined)))

## --- Per-cohort summary, z-scored WITHIN each cohort separately -------------
summarize_cohort <- function(data, ecm_genes) {
  populations <- sort(unique(data$cell_annotation))
  rows <- list()
  for (p in populations) {
    sub <- data[data$cell_annotation == p, ]
    n_cells <- nrow(sub)
    if (n_cells == 0) next
    for (g in ecm_genes) {
      vals <- sub[[g]]
      mean_expr <- mean(vals)
      pct_expr <- 100 * sum(vals > 0) / n_cells
      rows[[length(rows) + 1]] <- data.frame(population = p, gene = g, mean_expr = mean_expr,
                                              pct_expr = pct_expr, n_cells = n_cells)
    }
  }
  tbl <- do.call(rbind, rows)
  tbl$expr_z <- ave(tbl$mean_expr, tbl$gene, FUN = function(x) {
    if (length(x) < 2 || sd(x) == 0) return(rep(NA_real_, length(x)))
    (x - mean(x)) / sd(x)
  })
  tbl
}

all_rows <- list()
for (cohort in c("A", "B", "C")) {
  sub <- combined[combined$cohort == cohort, ]
  tbl <- summarize_cohort(sub, ecm_genes)
  tbl$cohort <- cohort
  all_rows[[cohort]] <- tbl
  cat(sprintf("\nCohort %s: %d populations with cells present\n", cohort, length(unique(tbl$population))))
  missing_pop <- setdiff(sort(unique(combined$cell_annotation)), unique(tbl$population))
  if (length(missing_pop) > 0) cat(sprintf("  No cells in: %s\n", paste(missing_pop, collapse = ", ")))
}
full_tbl <- do.call(rbind, all_rows)

tables_dir <- file.path(repo, "results", "tables")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
write.table(full_tbl, file.path(tables_dir, "ecm_panel_by_population_per_cohort.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

## --- 3-panel figure ---------------------------------------------------------
fig_dir <- file.path(repo, "results", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

full_tbl$gene <- factor(full_tbl$gene, levels = rev(ecm_genes))
full_tbl$cohort <- factor(full_tbl$cohort, levels = c("A", "B", "C"))
p <- ggplot(full_tbl[!is.na(full_tbl$mean_expr), ], aes(x = population, y = gene, size = pct_expr, color = expr_z)) +
  geom_point() +
  facet_wrap(~cohort, nrow = 1) +
  scale_color_gradient2(low = "blue", mid = "grey90", high = "red", midpoint = 0, name = "Relative\nexpression (z)\n(within cohort)") +
  scale_size_continuous(name = "% expressing", range = c(0, 8)) +
  labs(title = "ECM panel expression by population, split by cohort",
       x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(face = "bold", size = 12))
ggsave(file.path(fig_dir, "ecm_panel_by_population_per_cohort.png"), p, width = 15, height = 6, dpi = 150, bg = "white")

cat("\nWrote results/figures/ecm_panel_by_population_per_cohort.png\n")
cat("Done.\n")
