## Report figure: one combined QC overview across all 3 cohorts' UNFILTERED
## cells, 3 panels in one row, same style, points colored by QC OUTCOME
## (Kept / Removed: %MT / Removed: depth-complexity) using the exact
## per-cell criteria already applied in filter_cohort_a.R / filter_cohort_b.R
## / filter_cohort_c.R (via qc_utils.R::filter_seurat_cohort), PLUS each
## cohort's specimen-level minimum-cell drop, so the "Kept" fraction here
## exactly reproduces each cohort's true overall retention (A 87.4%, B 76.5%,
## C 75.2%) - this script only reads the existing seurat_{A,B,C}_raw.rds
## checkpoints' metadata and classifies/plots; it recomputes no QC and
## changes no threshold, filtering script, or downstream result.
##
## Thresholds replicated (NOT re-derived here - see the cited scripts):
##   A: nFeature_RNA > 500 & < 4500, nCount_RNA > 1000, percent_mt < 2
##   B: nFeature_RNA > 500 & < 4500, nCount_RNA > 1000, percent_mt < 10
##   C: nFeature_RNA > 100 & < 400, non-MT UMIs > 150, percent_mt < 95
##   all three: specimen dropped entirely if its POST-cell-filter count <= 20
##
## Outcome classification (3 categories, matching the task's priority rule -
## a cell failing both %MT and depth/complexity is classified %MT):
##   1. fails %MT threshold                       -> "Removed: %MT"
##   2. else fails depth/complexity threshold      -> "Removed: depth/complexity"
##   3. else cell itself passes, but its specimen
##      was dropped (<=20 cells post-filter)       -> "Removed: depth/complexity"
##      (folded in here: a specimen-level low-recovery drop is conceptually
##      a depth/complexity exclusion at the specimen level, not a new 4th
##      category - and it is empirically tiny: see the printed per-cohort
##      breakdown. C has 0 dropped specimens, so this never applies to C.)
##   4. else                                       -> "Kept"
##
## Input:  seurat_{A,B,C}_raw.rds (unfiltered, pre-QC checkpoints)
## Output: results/qc/qc_combined_ABC.png (7.2 x 2.4 in, 300 dpi)
##         results/qc/qc_combined_ABC_counts.tsv (per cohort x outcome, n and %)
##
## Run: Rscript R/qc_combined_ABC_figure.R

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
library(scales)

out_dir <- file.path(repo, "results", "qc")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

panel_labels <- c(A = "A (snRNA)", B = "B (scRNA)", C = "C (scRNA, derivative)")
outcome_levels <- c("Kept", "Removed: depth/complexity", "Removed: %MT")
outcome_colors <- c("Kept" = "grey80", "Removed: depth/complexity" = "#3366CC", "Removed: %MT" = "#E69F00")

MIN_SPECIMEN_CELLS <- 20

classify_cohort <- function(co) {
  cat(sprintf("\nCohort %s...\n", co))
  obj <- readRDS(file.path(SCRATCH_DIR, sprintf("seurat_%s_raw.rds", co)))
  m <- obj@meta.data
  rm(obj); gc(verbose = FALSE)

  df <- data.frame(cohort = co, sample_id = m$sample_id, nCount = m$nCount_RNA,
                    nFeature = m$nFeature_RNA, percent_mt = m$percent_mt)
  n_before <- nrow(df)
  df <- df[df$nCount > 0 & df$nFeature > 0, ]  # log-axis plotting hygiene only, see header
  n_dropped_log <- n_before - nrow(df)
  if (n_dropped_log > 0) cat(sprintf("  %d cells dropped (nCount/nFeature <= 0, unplottable on log axes)\n", n_dropped_log))

  if (co %in% c("A", "B")) {
    mt_max <- if (co == "A") 2 else 10
    df$fail_mt <- df$percent_mt >= mt_max
    df$fail_depth <- df$nCount <= 1000 | df$nFeature <= 500 | df$nFeature >= 4500
    df$cell_pass_qc <- !(df$fail_mt | df$fail_depth)
  } else {
    non_mt_umis <- df$nCount * (1 - df$percent_mt / 100)
    df$fail_mt <- df$percent_mt >= 95
    df$fail_depth <- df$nFeature <= 100 | df$nFeature >= 400 | non_mt_umis <= 150
    df$cell_pass_qc <- !(df$fail_mt | df$fail_depth)
  }

  post_filter_counts <- table(df$sample_id[df$cell_pass_qc])
  keep_specimens <- names(post_filter_counts)[post_filter_counts > MIN_SPECIMEN_CELLS]
  dropped_specimens <- setdiff(unique(df$sample_id), keep_specimens)
  df$specimen_dropped <- !(df$sample_id %in% keep_specimens)
  df$pass_final <- df$cell_pass_qc & !df$specimen_dropped
  cat(sprintf("  Specimens dropped (<=%d cells post-filter): %s\n", MIN_SPECIMEN_CELLS,
              if (length(dropped_specimens) > 0) paste(dropped_specimens, collapse = ", ") else "none"))

  n_specimen_drop_only <- sum(df$specimen_dropped & df$cell_pass_qc)
  cat(sprintf("  Cells passing per-cell QC but excluded via specimen-level drop: %d (%.2f%% of cohort) - folded into 'Removed: depth/complexity'\n",
              n_specimen_drop_only, 100 * n_specimen_drop_only / nrow(df)))

  df$outcome <- ifelse(df$pass_final, "Kept",
                 ifelse(df$fail_mt, "Removed: %MT", "Removed: depth/complexity"))
  df$outcome <- factor(df$outcome, levels = outcome_levels)

  n_kept <- sum(df$outcome == "Kept")
  cat(sprintf("  Overall retention: %d/%d = %.1f%%\n", n_kept, nrow(df), 100 * n_kept / nrow(df)))

  df
}

all_meta <- do.call(rbind, lapply(c("A", "B", "C"), classify_cohort))
all_meta$panel <- factor(panel_labels[all_meta$cohort], levels = panel_labels)

## --- Counts table: per cohort x outcome, n and % -----------------------
counts_tbl <- as.data.frame(table(cohort = all_meta$cohort, outcome = all_meta$outcome))
names(counts_tbl)[3] <- "n_cells"
totals <- aggregate(n_cells ~ cohort, data = counts_tbl, FUN = sum)
names(totals)[2] <- "n_total"
counts_tbl <- merge(counts_tbl, totals, by = "cohort")
counts_tbl$pct <- 100 * counts_tbl$n_cells / counts_tbl$n_total
counts_tbl <- counts_tbl[order(counts_tbl$cohort, counts_tbl$outcome), c("cohort", "outcome", "n_cells", "n_total", "pct")]
write.table(counts_tbl, file.path(out_dir, "qc_combined_ABC_counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== QC outcome counts (per cohort), for results/qc/qc_combined_ABC_counts.tsv ===\n")
print(counts_tbl, row.names = FALSE)
cat("\n=== Overall retention check (should match Table 1: A 87.4%, B 76.5%, C 75.2%) ===\n")
for (co in c("A", "B", "C")) {
  kept_pct <- counts_tbl$pct[counts_tbl$cohort == co & counts_tbl$outcome == "Kept"]
  cat(sprintf("  %s: %.1f%%\n", co, kept_pct))
}

## --- Threshold dashed lines, per panel (unchanged from the prior figure) ---
vlines <- data.frame(
  panel = factor(c(panel_labels["A"], panel_labels["B"]), levels = panel_labels),
  xintercept = c(1000, 1000)
)
hlines <- data.frame(
  panel = factor(c(panel_labels["A"], panel_labels["A"], panel_labels["B"], panel_labels["B"],
                    panel_labels["C"], panel_labels["C"]), levels = panel_labels),
  yintercept = c(500, 4500, 500, 4500, 100, 400)
)

## Plot order: randomly shuffled (not grouped by category), so no single
## outcome systematically draws on top of another at points of overlap -
## seed fixed and recorded for reproducibility.
SHUFFLE_SEED <- 42
set.seed(SHUFFLE_SEED)
cat(sprintf("\nShuffling plot order (seed = %d)\n", SHUFFLE_SEED))
plot_meta <- all_meta[sample(nrow(all_meta)), ]

## Per-panel QC-outcome percentage label (top-left), pulled directly from
## counts_tbl (the same numbers written to qc_combined_ABC_counts.tsv).
fmt_pct <- function(x) {
  ## whole-number style when the value is exactly an integer (e.g. "0%" for
  ## A's 0.0% depth/complexity), one decimal otherwise - matches the task's
  ## own example string.
  if (abs(x - round(x)) < 1e-9) sprintf("%d%%", round(x)) else sprintf("%.1f%%", x)
}
label_df <- data.frame(panel = factor(panel_labels, levels = panel_labels))
label_df$label <- vapply(names(panel_labels), function(co) {
  kept <- counts_tbl$pct[counts_tbl$cohort == co & counts_tbl$outcome == "Kept"]
  mt <- counts_tbl$pct[counts_tbl$cohort == co & counts_tbl$outcome == "Removed: %MT"]
  dc <- counts_tbl$pct[counts_tbl$cohort == co & counts_tbl$outcome == "Removed: depth/complexity"]
  sprintf("Kept %s\n%%MT %s | depth %s", fmt_pct(kept), fmt_pct(mt), fmt_pct(dc))
}, character(1))

p <- ggplot(plot_meta, aes(x = nCount, y = nFeature, color = outcome)) +
  geom_point(size = 0.3, alpha = 0.3, shape = 16) +
  geom_vline(data = vlines, aes(xintercept = xintercept), linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_hline(data = hlines, aes(yintercept = yintercept), linetype = "dashed", color = "black", linewidth = 0.3) +
  ## -Inf/Inf positions silently fail under scale_x_log10 (the transform is
  ## applied to -Inf before the panel-edge special-case kicks in, producing
  ## NaN and dropping the row) - use explicit finite data-space coordinates
  ## instead, valid since facet_wrap uses shared/fixed axis ranges across all
  ## 3 panels (same x/y window everywhere). Bottom-left, in the empty corner
  ## below the lowest dashed threshold line in every panel (C's is the
  ## highest floor at nFeature=100, still well above y=5) and left of the
  ## nCount=1000 vertical line in A/B (x=100, two short lines, comfortably
  ## inside one log decade).
  geom_text(data = label_df, aes(x = 100, y = 5, label = label), inherit.aes = FALSE,
            hjust = 0, vjust = 0, size = 1.7, color = "black", lineheight = 0.9) +
  scale_x_log10(breaks = 10^(2:5), labels = scales::label_log()) +
  scale_y_log10(breaks = 10^(1:4), labels = scales::label_log()) +
  scale_color_manual(values = outcome_colors, name = "QC outcome", breaks = outcome_levels) +
  guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1))) +
  facet_wrap(~panel, nrow = 1) +
  labs(x = "nCount (log10)", y = "nFeature (log10)") +
  theme_minimal(base_size = 7) +
  theme(strip.text = element_text(face = "bold", size = 7),
        legend.key.height = unit(0.3, "cm"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 5.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.spacing.x = unit(0.9, "lines"))

out_path <- file.path(out_dir, "qc_combined_ABC.png")
ggsave(out_path, p, width = 7.2, height = 2.4, dpi = 300, bg = "white")
cat(sprintf("\nWrote %s\n", out_path))
cat("Done.\n")
