## Cohort A QC: per-nucleus metrics computed on the author-annotated set
## (see load_cohort_a.R) - NOT a fresh droplet-style filter. DATA_DICTIONARY
## is explicit that droplet-cell thresholds must not be applied uniformly to
## nuclei, and that the author's analyzed set already IS the inclusion
## criterion here. This script characterizes that set (and flags any
## specimen with conspicuously low recovery) rather than re-filtering it.
##
## Run: Rscript R/qc_cohort_a.R  (needs cohort_a_loaded.rds - see
## R/run_load_cohort_a.R)

.args <- commandArgs(trailingOnly = FALSE)
.script_path <- sub("--file=", "", .args[grep("--file=", .args)])
repo <- normalizePath(file.path(dirname(.script_path), ".."))
readRenviron(file.path(repo, ".Renviron"))
source(file.path(repo, "config.R"))
source(file.path(repo, "R", "qc_utils.R"))
library(ggplot2)

sample_manifest <- read.delim(SAMPLE_MANIFEST_PATH, stringsAsFactors = FALSE)
a_manifest <- sample_manifest[sample_manifest$cohort == "A", , drop = FALSE]

cache_path <- file.path(SCRATCH_DIR, "cohort_a_loaded.rds")
if (!file.exists(cache_path)) {
  stop("cohort_a_loaded.rds not found - run R/run_load_cohort_a.R first.")
}
cohort_a <- readRDS(cache_path)

cat("Computing per-nucleus QC metrics for", length(cohort_a), "specimens...\n")
qc_list <- vector("list", length(cohort_a))
for (sid in names(cohort_a)) {
  r <- cohort_a[[sid]]
  mt_ids <- mt_row_ids(r$features, id_col = "gene_id")
  m <- compute_qc_metrics(r$matrix, mt_ids)
  m$sample_id <- sid
  m$broad_label <- r$cell_metadata$Cell_type.shorter
  m$detailed_label <- r$cell_metadata$Cell_type.detailed
  qc_list[[sid]] <- m
}
qc <- do.call(rbind, qc_list)
qc <- merge(qc, a_manifest[, c("sample_id", "donor_id", "condition")], by = "sample_id")

cat("Total nuclei:", nrow(qc), "\n")
cat("MT genes identified (by symbol, cohort A features):",
    length(mt_row_ids(cohort_a[[1]]$features, "gene_id")), "\n")

## Cache the per-cell QC table itself (cheap, ~5MB), not just the plots/
## summary - so a figure-only change later doesn't require re-running the
## ~20 min full matrix load again.
saveRDS(qc, file.path(SCRATCH_DIR, "qc_cohort_a.rds"))

fig_dir <- file.path(repo, "results", "figures", "A")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## --- Per-specimen summary table, flagging low-recovery specimens ---------
per_specimen <- aggregate(
  cbind(n_count, n_feature, percent_mt) ~ sample_id + donor_id + condition,
  data = qc, FUN = median
)
n_cells_tab <- as.data.frame(table(sample_id = qc$sample_id))
colnames(n_cells_tab) <- c("sample_id", "n_cells")
per_specimen <- merge(per_specimen, n_cells_tab, by = "sample_id")
colnames(per_specimen)[colnames(per_specimen) %in% c("n_count", "n_feature", "percent_mt")] <-
  c("median_n_count", "median_n_feature", "median_percent_mt")
per_specimen <- per_specimen[order(per_specimen$n_cells), ]
# Flagged, not excluded - a low-recovery flag informs downstream population-
# choice/robustness decisions rather than dropping the specimen outright.
per_specimen$low_recovery_flag <- per_specimen$n_cells < LOW_RECOVERY_THRESHOLD

write.table(
  per_specimen, file.path(repo, "results", "tables", "qc_cohort_a_per_specimen.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
cat(sprintf("\nSpecimens flagged for low nucleus recovery (<%d):\n", LOW_RECOVERY_THRESHOLD))
print(per_specimen[per_specimen$low_recovery_flag, c("sample_id", "donor_id", "condition", "n_cells")],
      row.names = FALSE)

## --- Figures ---------------------------------------------------------------
qc$condition <- factor(qc$condition, levels = c("Tumor", "Adjacent"))
specimen_order <- per_specimen$sample_id[order(per_specimen$condition, per_specimen$n_cells)]
qc$sample_id <- factor(qc$sample_id, levels = specimen_order)

p1 <- ggplot(qc, aes(x = sample_id, y = n_count, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE) +
  scale_y_log10() +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(title = "Cohort A: UMI counts per nucleus, by specimen", x = NULL, y = "n_count (log10)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))
ggsave(file.path(fig_dir, "qc_cohort_a_ncount_violin.png"), p1, width = 12, height = 5, dpi = 150)

p2 <- ggplot(qc, aes(x = sample_id, y = percent_mt, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE) +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(title = "Cohort A: % mitochondrial counts per nucleus, by specimen", x = NULL, y = "% MT") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))
ggsave(file.path(fig_dir, "qc_cohort_a_percentmt_violin.png"), p2, width = 12, height = 5, dpi = 150)

p2b <- ggplot(qc, aes(x = sample_id, y = n_feature, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE) +
  scale_y_log10() +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(title = "Cohort A: genes detected per nucleus, by specimen", x = NULL, y = "n_feature (log10)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))
ggsave(file.path(fig_dir, "qc_cohort_a_nfeature_violin.png"), p2b, width = 12, height = 5, dpi = 150)

p3 <- ggplot(qc, aes(x = n_count, y = n_feature, color = condition)) +
  geom_point(alpha = 0.15, size = 0.4) +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = CONDITION_COLORS) +
  labs(title = "Cohort A: n_count vs n_feature per nucleus", x = "n_count (log10)", y = "n_feature (log10)") +
  theme_minimal()
ggsave(file.path(fig_dir, "qc_cohort_a_count_vs_feature.png"), p3, width = 7, height = 6, dpi = 150)

p3b <- linear_count_vs_feature_plot(qc, "Cohort A: n_count vs n_feature per nucleus (linear axes)")
ggsave(file.path(fig_dir, "qc_cohort_a_count_vs_feature_linear.png"), p3b, width = 7, height = 6, dpi = 150)

per_specimen$condition <- factor(per_specimen$condition, levels = c("Tumor", "Adjacent"))
p4 <- ggplot(per_specimen, aes(x = reorder(sample_id, n_cells), y = n_cells, fill = condition)) +
  geom_col() +
  scale_fill_manual(values = CONDITION_COLORS) +
  coord_flip() +
  labs(title = "Cohort A: nuclei recovered per specimen", x = NULL, y = "n_cells") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6))
ggsave(file.path(fig_dir, "qc_cohort_a_cells_per_specimen.png"), p4, width = 7, height = 8, dpi = 150)

cat("\nWrote 6 figures and 1 table to results/figures and results/tables.\n")
