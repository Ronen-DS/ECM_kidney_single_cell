## Verifies the final cell annotation (merged_ABC_annotated.rds) using
## canonical lineage markers NOT in the fixed ECM panel (COL1A1, COL1A2,
## COL3A1, COL5A1, COL6A1, COL6A2, FN1, SPARC, DCN, LUM). The integrated
## object only carries the 2000 SelectIntegrationFeatures genes (chosen by
## cross-cohort variability, not curated markers), so marker expression is
## pulled from each cohort's own FULL normalized checkpoint instead, matched
## back to the final annotation/cluster by cell id.
##
## A marker absent from a cohort's own gene reference (e.g. PECAM1 in B) is
## excluded from that cohort's cells for that marker specifically (genuine
## NA, not a zero) - the per-label/per-cluster averages are computed only
## over cells from cohorts where the gene was actually measured.
##
## "Relative expression" is reported as each marker's per-label (or
## per-cluster) mean expression scaled (z-scored) across labels/clusters -
## NOT a statistical enrichment test (no hypergeometric/GSEA-style
## over-representation), just a z-scored mean - alongside % of cells
## expressing, the standard dot-plot pairing, chosen specifically because
## single-cell marker dropout (especially in cohort C) makes raw per-cell
## detection unreliable on its own.
##
## Input:  merged_ABC_annotated.rds, seurat_{A,B,C}_normalized.rds
## Output: results/tables/integration/marker_verification_by_label.tsv
##         results/tables/integration/marker_verification_by_cluster.tsv
##         results/figures/integration/marker_dotplot_by_label.png
##         results/figures/integration/marker_dotplot_by_cluster.png
##
## Run: Rscript R/verify_annotation_markers.R

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

markers <- c(
  "PECAM1", "CDH5",       # Endothelial
  "CD68", "LYZ",          # Myeloid
  "CD3D", "PTPRC",        # Lymphoid (PTPRC = pan-immune cross-check)
  "CA9", "NDUFA4L2",      # Malignant_epithelial (ccRCC-specific)
  "EPCAM",                # shared epithelial lineage
  "RGS5", "PDGFRB",       # Pericyte
  "PDGFRA", "ACTA2"       # Stromal_fibro_myo. TAGLN was also tried as a
  # sharper myofibroblast-specific alternative to ACTA2, but it landed highest
  # in Pericyte too (same as ACTA2) - confirming the Pericyte/Stromal_fibro_myo
  # overlap is real shared contractile-lineage biology here, not an artifact
  # of this particular marker choice. Dropped to avoid redundant, non-resolving
  # entries in the panel.
)

cat("Loading merged_ABC_annotated.rds for metadata only...\n")
obj_meta <- readRDS(file.path(SCRATCH_DIR, "merged_ABC_annotated.rds"))
meta <- data.frame(
  cell_key = colnames(obj_meta),
  cohort = obj_meta$cohort,
  cell_annotation = obj_meta$cell_annotation,
  cluster = obj_meta$rpca_snn_res.1,
  stringsAsFactors = FALSE
)
rm(obj_meta)
gc(verbose = FALSE)
cat(sprintf("  %d cells\n", nrow(meta)))

cat("\n=== Marker presence per cohort ===\n")
marker_frames <- list()
for (cohort in c("A", "B", "C")) {
  cat(sprintf("\nCohort %s...\n", cohort))
  obj <- readRDS(file.path(SCRATCH_DIR, sprintf("seurat_%s_normalized.rds", cohort)))
  present <- intersect(markers, rownames(obj))
  absent <- setdiff(markers, rownames(obj))
  cat(sprintf("  Present: %s\n", paste(present, collapse = ", ")))
  if (length(absent) > 0) cat(sprintf("  ABSENT (excluded for this cohort, not zero-filled): %s\n", paste(absent, collapse = ", ")))

  mat <- GetAssayData(obj, assay = "RNA", layer = "data")[present, , drop = FALSE]
  df <- as.data.frame(as.matrix(t(mat)))
  df$cell_key <- colnames(obj)
  for (m in absent) df[[m]] <- NA_real_ # explicit, genuine NA - never measured in this cohort
  marker_frames[[cohort]] <- df[, c("cell_key", markers)]
  rm(obj, mat)
  gc(verbose = FALSE)
}
marker_df <- do.call(rbind, marker_frames)
cat(sprintf("\nCombined marker table: %d cells x %d markers\n", nrow(marker_df), length(markers)))

combined <- merge(meta, marker_df, by = "cell_key")
cat(sprintf("After matching to annotation: %d cells\n", nrow(combined)))

## --- Summarize: mean expression (NA-excluded) + %expressing, by group -------
summarize_by <- function(data, group_col, markers) {
  groups <- sort(unique(data[[group_col]]))
  rows <- list()
  for (g in groups) {
    sub <- data[data[[group_col]] == g, ]
    for (m in markers) {
      vals <- sub[[m]]
      n_measured <- sum(!is.na(vals))
      mean_expr <- if (n_measured > 0) mean(vals, na.rm = TRUE) else NA_real_
      pct_expr <- if (n_measured > 0) 100 * sum(vals > 0, na.rm = TRUE) / n_measured else NA_real_
      rows[[length(rows) + 1]] <- data.frame(
        group = g, marker = m, mean_expr = mean_expr, pct_expr = pct_expr, n_measured = n_measured
      )
    }
  }
  tbl <- do.call(rbind, rows)
  ## relative expression: z-score each marker's mean_expr across groups
  ## (a standardized mean, not a statistical enrichment test)
  tbl$expr_z <- ave(tbl$mean_expr, tbl$marker, FUN = function(x) {
    if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) return(rep(NA_real_, length(x)))
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  })
  tbl
}

cat("\nSummarizing by assigned label...\n")
by_label <- summarize_by(combined, "cell_annotation", markers)
cat("\nSummarizing by cluster (resolution 1)...\n")
combined$cluster <- as.character(combined$cluster)
by_cluster <- summarize_by(combined, "cluster", markers)

tables_dir <- file.path(repo, "results", "tables", "integration")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
write.table(by_label, file.path(tables_dir, "marker_verification_by_label.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(by_cluster, file.path(tables_dir, "marker_verification_by_cluster.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== By label: dominant (highest relative expression z) label per marker ===\n")
for (m in markers) {
  sub <- by_label[by_label$marker == m & !is.na(by_label$expr_z), ]
  if (nrow(sub) == 0) { cat(sprintf("%-10s: no data\n", m)); next }
  top <- sub[which.max(sub$expr_z), ]
  cat(sprintf("%-10s: highest in '%s' (mean=%.3f, z=%.2f, %%expr=%.1f%%, n=%d)\n",
              m, top$group, top$mean_expr, top$expr_z, top$pct_expr, top$n_measured))
}

## --- Dot plots ------------------------------------------------------------
fig_dir <- file.path(repo, "results", "figures", "integration")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

by_label$marker <- factor(by_label$marker, levels = rev(markers))
p1 <- ggplot(by_label[!is.na(by_label$mean_expr), ], aes(x = group, y = marker, size = pct_expr, color = expr_z)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "grey90", high = "red", midpoint = 0, name = "Relative\nexpression (z)") +
  scale_size_continuous(name = "% expressing", range = c(0, 8)) +
  labs(title = "Canonical marker verification by assigned cell population",
       x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, "marker_dotplot_by_label.png"), p1, width = 9, height = 6, dpi = 150, bg = "white")

by_cluster$cluster_num <- as.numeric(by_cluster$group)
by_cluster$marker <- factor(by_cluster$marker, levels = rev(markers))
p2 <- ggplot(by_cluster[!is.na(by_cluster$mean_expr), ], aes(x = reorder(group, cluster_num), y = marker, size = pct_expr, color = expr_z)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "grey90", high = "red", midpoint = 0, name = "Relative\nexpression (z)") +
  scale_size_continuous(name = "% expressing", range = c(0, 6)) +
  labs(title = "Canonical marker verification by resolution-1 cluster",
       x = "cluster", y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7))
ggsave(file.path(fig_dir, "marker_dotplot_by_cluster.png"), p2, width = 14, height = 6, dpi = 150, bg = "white")

cat("\nWrote marker_dotplot_by_label.png and marker_dotplot_by_cluster.png\n")
cat("Done.\n")
