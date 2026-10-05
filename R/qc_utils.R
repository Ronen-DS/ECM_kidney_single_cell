## Shared QC helpers used by all three per-cohort QC scripts.
##
## Per-cell metrics: nCount (total UMI/reads), nFeature (genes detected),
## percent_mt (% of counts from mitochondrial genes, standard "MT-" symbol
## prefix). Mitochondrial genes are identified by SYMBOL (the only feature
## description every cohort shares), then mapped back to whatever row
## identifier that cohort's matrix actually uses (gene_id for A/C, symbol
## itself for B) - this keeps the same logic correct for all three despite
## their different row-identifier conventions.

library(Matrix)

## Fixed, explicit Tumor/Adjacent color mapping shared by every figure in
## every QC script (matches ggplot's default 2-color hue palette). Used via
## scale_fill/color_manual(values = CONDITION_COLORS) rather than relying on
## factor level order, which silently flips between a violin plot (built
## from a releveled factor) and a bar chart (built from aggregate() output
## that was never releveled) - exactly the inconsistency this fixes.
CONDITION_COLORS <- c(Tumor = "#F8766D", Adjacent = "#00BFC4")

#' Row identifiers (whatever they are - gene_id or symbol) corresponding to
#' mitochondrial genes, given a features table with a gene_symbol column and
#' a column of the identifier actually used as matrix rownames.
mt_row_ids <- function(features_df, id_col, symbol_col = "gene_symbol", pattern = "^MT-") {
  is_mt <- grepl(pattern, features_df[[symbol_col]])
  unique(features_df[[id_col]][is_mt])
}

#' Symbols that map to more than one distinct gene_id within one cohort's
#' features table - see the loader duplicate-symbol notes. Any symbol-keyed
#' cross-cohort lookup should exclude these rather than guess which row is
#' "the" gene.
ambiguous_symbols <- function(features_df, symbol_col = "gene_symbol") {
  syms <- features_df[[symbol_col]]
  unique(syms[duplicated(syms)])
}

## Low-recovery flag threshold, shared across all 3 cohorts' per-specimen
## summary tables and "cells per specimen" bar charts (dashed reference
## line). A flag, not an exclusion - see each qc_cohort_*.R script.
LOW_RECOVERY_THRESHOLD <- 1000

#' A linear-axis (not log10) n_count vs n_feature scatter, for readers who
#' want to see the actual scale rather than a log-compressed one. Clips the
#' visible range to the given percentile (default 99th) of each axis via
#' coord_cartesian - this only changes what's IN FRAME, it does not drop or
#' alter any data point, so a few extreme cells don't compress the entire
#' rest of the distribution into a corner of the plot.
linear_count_vs_feature_plot <- function(qc, title, clip_quantile = 0.99) {
  xlim <- c(0, stats::quantile(qc$n_count, clip_quantile, na.rm = TRUE))
  ylim <- c(0, stats::quantile(qc$n_feature, clip_quantile, na.rm = TRUE))
  ggplot2::ggplot(qc, ggplot2::aes(x = n_count, y = n_feature, color = condition)) +
    ggplot2::geom_point(alpha = 0.15, size = 0.4) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
    ggplot2::scale_color_manual(values = CONDITION_COLORS) +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("View clipped at the %.0fth percentile; no points removed", 100 * clip_quantile),
      x = "n_count", y = "n_feature"
    ) +
    ggplot2::theme_minimal()
}

#' Per-cell nFeature_RNA/nCount_RNA/percent_mt filter, then a specimen-level
#' minimum-cell filter, applied to a Seurat object's own metadata columns.
#' min_specimen_cells is checked against each specimen's RESULTING
#' (post-cell-filter) cell count, not its raw count - a specimen that has
#' plenty of raw cells but few surviving the per-cell filters is dropped
#' entirely (its surviving handful, if any, is discarded too), rather than
#' kept with a thin remainder. Returns the filtered object plus a
#' before/after summary.
#'
#' non_mt_umis_min is an OPTIONAL extra criterion (default NULL = unused, so
#' existing callers are unaffected): non-MT UMIs (nCount * (1 - percent_mt/100),
#' the actual signal left once mitochondrial reads are set aside) must exceed
#' this value. It mathematically subsumes an nCount floor of the same value
#' (non-MT UMIs <= nCount always), so a cohort using this criterion typically
#' passes n_count_min = -Inf rather than double-applying the same idea -
#' introduced for cohort C, where %MT and raw depth turned out to be only
#' loosely coupled (see derive_c_mt_threshold.R).
filter_seurat_cohort <- function(obj, min_specimen_cells,
                                  n_feature_min, n_feature_max,
                                  n_count_min, n_count_max,
                                  percent_mt_max,
                                  non_mt_umis_min = NULL,
                                  sample_col = "sample_id") {
  meta <- obj@meta.data
  n_specimens_before <- length(unique(meta[[sample_col]]))

  cell_pass_qc <- meta$nFeature_RNA > n_feature_min & meta$nFeature_RNA < n_feature_max &
    meta$nCount_RNA > n_count_min & meta$nCount_RNA < n_count_max &
    meta$percent_mt < percent_mt_max

  if (!is.null(non_mt_umis_min)) {
    non_mt_umis <- meta$nCount_RNA * (1 - meta$percent_mt / 100)
    cell_pass_qc <- cell_pass_qc & non_mt_umis > non_mt_umis_min
  }

  post_filter_counts <- table(meta[[sample_col]][cell_pass_qc])
  keep_specimens <- names(post_filter_counts)[post_filter_counts > min_specimen_cells]
  dropped_specimens <- setdiff(unique(meta[[sample_col]]), keep_specimens)

  cell_pass <- cell_pass_qc & meta[[sample_col]] %in% keep_specimens

  list(
    obj = obj[, cell_pass],
    dropped_specimens = dropped_specimens,
    n_cells_before = nrow(meta),
    n_cells_after = sum(cell_pass),
    n_specimens_before = n_specimens_before,
    n_specimens_after = length(keep_specimens)
  )
}

#' Writes the same 6 standard QC figures used by the raw qc_cohort_*.R
#' scripts (3 violins, log-log scatter, linear scatter, cells-per-specimen
#' bar), given a qc data frame (needs sample_id/condition/n_count/n_feature/
#' percent_mt) and its per-specimen summary (needs sample_id/condition/
#' n_cells). Sharing this between raw and post-filter QC keeps the two sets
#' of figures directly comparable (same style, same file names, different
#' directory).
write_standard_qc_figures <- function(qc, per_specimen, fig_dir, cohort_label, file_prefix,
                                       threshold_line = LOW_RECOVERY_THRESHOLD,
                                       violin_width = 9, violin_height = 5, bar_height = 6) {
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  qc$condition <- factor(qc$condition, levels = c("Tumor", "Adjacent"))
  specimen_order <- per_specimen$sample_id[order(per_specimen$condition, per_specimen$n_cells)]
  qc$sample_id <- factor(qc$sample_id, levels = specimen_order)
  per_specimen$condition <- factor(per_specimen$condition, levels = c("Tumor", "Adjacent"))

  p1 <- ggplot2::ggplot(qc, ggplot2::aes(x = sample_id, y = n_count, fill = condition)) +
    ggplot2::geom_violin(scale = "width", trim = TRUE) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_fill_manual(values = CONDITION_COLORS) +
    ggplot2::labs(title = sprintf("%s: UMI counts per cell, by specimen", cohort_label), x = NULL, y = "n_count (log10)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))
  ggplot2::ggsave(file.path(fig_dir, paste0(file_prefix, "_ncount_violin.png")), p1, width = violin_width, height = violin_height, dpi = 150)

  p2 <- ggplot2::ggplot(qc, ggplot2::aes(x = sample_id, y = n_feature, fill = condition)) +
    ggplot2::geom_violin(scale = "width", trim = TRUE) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_fill_manual(values = CONDITION_COLORS) +
    ggplot2::labs(title = sprintf("%s: genes detected per cell, by specimen", cohort_label), x = NULL, y = "n_feature (log10)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))
  ggplot2::ggsave(file.path(fig_dir, paste0(file_prefix, "_nfeature_violin.png")), p2, width = violin_width, height = violin_height, dpi = 150)

  p3 <- ggplot2::ggplot(qc, ggplot2::aes(x = sample_id, y = percent_mt, fill = condition)) +
    ggplot2::geom_violin(scale = "width", trim = TRUE) +
    ggplot2::scale_fill_manual(values = CONDITION_COLORS) +
    ggplot2::labs(title = sprintf("%s: %% mitochondrial counts per cell, by specimen", cohort_label), x = NULL, y = "% MT") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))
  ggplot2::ggsave(file.path(fig_dir, paste0(file_prefix, "_percentmt_violin.png")), p3, width = violin_width, height = violin_height, dpi = 150)

  p4 <- ggplot2::ggplot(qc, ggplot2::aes(x = n_count, y = n_feature, color = condition)) +
    ggplot2::geom_point(alpha = 0.15, size = 0.4) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::scale_color_manual(values = CONDITION_COLORS) +
    ggplot2::labs(title = sprintf("%s: n_count vs n_feature per cell", cohort_label), x = "n_count (log10)", y = "n_feature (log10)") +
    ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(fig_dir, paste0(file_prefix, "_count_vs_feature.png")), p4, width = 7, height = 6, dpi = 150)

  p5 <- linear_count_vs_feature_plot(qc, sprintf("%s: n_count vs n_feature per cell (linear axes)", cohort_label))
  ggplot2::ggsave(file.path(fig_dir, paste0(file_prefix, "_count_vs_feature_linear.png")), p5, width = 7, height = 6, dpi = 150)

  p6 <- ggplot2::ggplot(per_specimen, ggplot2::aes(x = reorder(sample_id, n_cells), y = n_cells, fill = condition)) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = threshold_line, linetype = "dashed", color = "red") +
    ggplot2::scale_fill_manual(values = CONDITION_COLORS) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = sprintf("%s: cells per specimen", cohort_label), x = NULL, y = "n_cells") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6))
  ggplot2::ggsave(file.path(fig_dir, paste0(file_prefix, "_cells_per_specimen.png")), p6, width = 7, height = bar_height, dpi = 150)

  invisible(NULL)
}

#' Computes the 3 standard per-cell QC metrics for one specimen's matrix
#' (features x cells, any sparse Matrix). mt_ids must already be in the same
#' identifier space as rownames(mat) (see mt_row_ids()).
compute_qc_metrics <- function(mat, mt_ids) {
  n_count <- Matrix::colSums(mat)
  n_feature <- Matrix::colSums(mat > 0)
  mt_ids_present <- intersect(mt_ids, rownames(mat))
  mt_count <- if (length(mt_ids_present) > 0) {
    Matrix::colSums(mat[mt_ids_present, , drop = FALSE])
  } else {
    rep(0, ncol(mat))
  }
  percent_mt <- ifelse(n_count > 0, 100 * mt_count / n_count, NA_real_)
  data.frame(
    cell_key = colnames(mat),
    n_count = n_count,
    n_feature = n_feature,
    percent_mt = percent_mt,
    row.names = NULL
  )
}
