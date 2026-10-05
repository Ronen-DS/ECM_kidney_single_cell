## Compact overview of ECM-panel expression across the major cell
## populations (Section 2). Mirrors verify_annotation_markers.R's approach:
## expression is pulled from each cohort's own FULL normalized checkpoint
## (the integrated object only carries the 2000 SelectIntegrationFeatures
## genes, not necessarily these 10), matched back to the final cell
## annotation (merged_ABC_annotated.rds) by cell id. Presence of all 10 ECM
## panel genes in every cohort is verified explicitly here rather than
## assumed from the earlier DATA_DICTIONARY-based check.
##
## "Relative expression" is each gene's per-population mean log-normalized
## expression z-scored across populations (not a statistical enrichment
## test) - the same approach and terminology used for the lineage-marker
## verification, for consistency.
##
## Input:  merged_ABC_annotated.rds, seurat_{A,B,C}_normalized.rds
## Output: results/tables/ecm_panel_by_population.tsv
##         results/figures/ecm_panel_by_population.png
##
## Run: Rscript R/ecm_panel_by_population.R

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

cat("\n=== ECM panel gene presence per cohort ===\n")
ecm_frames <- list()
for (cohort in c("A", "B", "C")) {
  cat(sprintf("\nCohort %s...\n", cohort))
  obj <- readRDS(file.path(SCRATCH_DIR, sprintf("seurat_%s_normalized.rds", cohort)))
  present <- intersect(ecm_genes, rownames(obj))
  absent <- setdiff(ecm_genes, rownames(obj))
  cat(sprintf("  Present: %d/%d (%s)\n", length(present), length(ecm_genes), paste(present, collapse = ", ")))
  if (length(absent) > 0) cat(sprintf("  ABSENT (excluded for this cohort, not zero-filled): %s\n", paste(absent, collapse = ", ")))

  mat <- GetAssayData(obj, assay = "RNA", layer = "data")[present, , drop = FALSE]
  df <- as.data.frame(as.matrix(t(mat)))
  df$cell_key <- colnames(obj)
  for (g in absent) df[[g]] <- NA_real_ # explicit, genuine NA - never measured in this cohort
  ecm_frames[[cohort]] <- df[, c("cell_key", ecm_genes)]
  rm(obj, mat)
  gc(verbose = FALSE)
}
ecm_df <- do.call(rbind, ecm_frames)
cat(sprintf("\nCombined ECM table: %d cells x %d genes\n", nrow(ecm_df), length(ecm_genes)))

combined <- merge(meta, ecm_df, by = "cell_key")
cat(sprintf("After matching to annotation: %d cells\n", nrow(combined)))

## --- Summarize: mean expression (NA-excluded) + %expressing, by population ---
populations <- sort(unique(combined$cell_annotation))
rows <- list()
for (p in populations) {
  sub <- combined[combined$cell_annotation == p, ]
  for (g in ecm_genes) {
    vals <- sub[[g]]
    n_measured <- sum(!is.na(vals))
    mean_expr <- if (n_measured > 0) mean(vals, na.rm = TRUE) else NA_real_
    pct_expr <- if (n_measured > 0) 100 * sum(vals > 0, na.rm = TRUE) / n_measured else NA_real_
    rows[[length(rows) + 1]] <- data.frame(
      population = p, gene = g, mean_expr = mean_expr, pct_expr = pct_expr,
      n_cells = nrow(sub), n_measured = n_measured
    )
  }
}
tbl <- do.call(rbind, rows)
## relative expression: z-score each gene's mean_expr across populations
## (a standardized mean, not a statistical enrichment test)
tbl$expr_z <- ave(tbl$mean_expr, tbl$gene, FUN = function(x) {
  if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
})

tables_dir <- file.path(repo, "results", "tables")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
write.table(tbl, file.path(tables_dir, "ecm_panel_by_population.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== Per-gene summary: population with highest relative expression ===\n")
for (g in ecm_genes) {
  sub <- tbl[tbl$gene == g & !is.na(tbl$expr_z), ]
  if (nrow(sub) == 0) { cat(sprintf("%-10s: no data\n", g)); next }
  top <- sub[which.max(sub$expr_z), ]
  cat(sprintf("%-10s: highest in '%s' (mean=%.3f, z=%.2f, %%expr=%.1f%%, n=%d)\n",
              g, top$population, top$mean_expr, top$expr_z, top$pct_expr, top$n_measured))
}

## --- Compact dot plot -------------------------------------------------------
fig_dir <- file.path(repo, "results", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

tbl$gene <- factor(tbl$gene, levels = rev(ecm_genes))
p <- ggplot(tbl[!is.na(tbl$mean_expr), ], aes(x = population, y = gene, size = pct_expr, color = expr_z)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "grey90", high = "red", midpoint = 0, name = "Relative\nexpression (z)") +
  scale_size_continuous(name = "% expressing", range = c(0, 10)) +
  labs(title = "ECM panel expression across major cell populations (A+B+C)",
       x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, "ecm_panel_by_population.png"), p, width = 9, height = 6, dpi = 150, bg = "white")

cat("\nWrote results/figures/ecm_panel_by_population.png\n")
cat("Done.\n")
