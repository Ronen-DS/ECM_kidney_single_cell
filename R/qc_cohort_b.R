## Cohort B QC: reads the cached Seurat checkpoint (see build_seurat_b.R -
## the actual load-26-CSVs step, ~35 min, lives there now, not here) and
## produces the same 4 standard QC figures as cohorts A/C (n_count violin,
## n_feature violin, %MT violin, n_count vs n_feature scatter), plus B's
## own extras: a check that our recomputed nCount/nFeature/%MT are
## consistent with the SUPPLIED author QC fields (TotalUMI, gene.num,
## MT.ratio) - per DATA_DICTIONARY, these are fields to validate, not a
## license to derive a fresh rejection rule - and the supplied doublet score.
##
## Run: Rscript R/qc_cohort_b.R  (needs seurat_B_raw.rds - see R/build_seurat_b.R)

.args <- commandArgs(trailingOnly = FALSE)
.script_path <- sub("--file=", "", .args[grep("--file=", .args)])
repo <- normalizePath(file.path(dirname(.script_path), ".."))
readRenviron(file.path(repo, ".Renviron"))
source(file.path(repo, "config.R"))
source(file.path(repo, "R", "qc_utils.R"))
library(Seurat)
library(ggplot2)

cache_path <- file.path(SCRATCH_DIR, "seurat_B_raw.rds")
if (!file.exists(cache_path)) {
  stop("seurat_B_raw.rds not found - run R/build_seurat_b.R first.")
}
cat("Loading cached Seurat object...\n")
obj <- readRDS(cache_path)

qc <- obj@meta.data
qc$n_count <- qc$nCount_RNA
qc$n_feature <- qc$nFeature_RNA
cat("Total cells:", nrow(qc), "\n")

## --- Consistency check: recomputed vs. supplied QC fields ------------------
cor_count <- cor(qc$n_count, qc$TotalUMI, use = "complete.obs")
cor_feature <- cor(qc$n_feature, qc$gene.num, use = "complete.obs")
cor_mt <- cor(qc$percent_mt, qc$MT.ratio_pct, use = "complete.obs")
cat(sprintf("\nCorrelation (recomputed vs supplied): n_count~TotalUMI r=%.4f | n_feature~gene.num r=%.4f | percent_mt~MT.ratio r=%.4f\n",
            cor_count, cor_feature, cor_mt))

consistency_summary <- data.frame(
  metric = c("n_count vs TotalUMI", "n_feature vs gene.num", "percent_mt vs MT.ratio(*100)"),
  pearson_r = c(cor_count, cor_feature, cor_mt),
  median_abs_diff = c(
    median(abs(qc$n_count - qc$TotalUMI), na.rm = TRUE),
    median(abs(qc$n_feature - qc$gene.num), na.rm = TRUE),
    median(abs(qc$percent_mt - qc$MT.ratio_pct), na.rm = TRUE)
  )
)
write.table(consistency_summary, file.path(repo, "results", "tables", "qc_cohort_b_consistency.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
print(consistency_summary, row.names = FALSE)

## --- Per-specimen summary table ---------------------------------------------
per_specimen <- aggregate(
  cbind(n_count, n_feature, percent_mt, doublest.score) ~ sample_id + donor_id + condition,
  data = qc, FUN = median
)
n_cells_tab <- as.data.frame(table(sample_id = qc$sample_id))
colnames(n_cells_tab) <- c("sample_id", "n_cells")
per_specimen <- merge(per_specimen, n_cells_tab, by = "sample_id")
colnames(per_specimen)[colnames(per_specimen) %in% c("n_count", "n_feature", "percent_mt", "doublest.score")] <-
  c("median_n_count", "median_n_feature", "median_percent_mt", "median_doublet_score")
per_specimen <- per_specimen[order(per_specimen$n_cells), ]
per_specimen$low_recovery_flag <- per_specimen$n_cells < LOW_RECOVERY_THRESHOLD
write.table(per_specimen, file.path(repo, "results", "tables", "qc_cohort_b_per_specimen.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nn cells with missing/unresolved ano.l1:", sum(is.na(qc$ano.l1)), "\n")

## --- Figures -----------------------------------------------------------------
fig_dir <- file.path(repo, "results", "figures", "B")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

qc$condition <- factor(qc$condition, levels = c("Tumor", "Adjacent"))
specimen_order <- per_specimen$sample_id[order(per_specimen$condition, per_specimen$n_cells)]
qc$sample_id <- factor(qc$sample_id, levels = specimen_order)
per_specimen$condition <- factor(per_specimen$condition, levels = c("Tumor", "Adjacent"))

## The 4 standard figures, matching cohorts A and C exactly -------------------
p1 <- ggplot(qc, aes(x = sample_id, y = n_count, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE) +
  scale_y_log10() +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(title = "Cohort B: UMI counts per cell, by specimen", x = NULL, y = "n_count (log10)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))
ggsave(file.path(fig_dir, "qc_cohort_b_ncount_violin.png"), p1, width = 10, height = 5, dpi = 150)

p1b <- ggplot(qc, aes(x = sample_id, y = n_feature, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE) +
  scale_y_log10() +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(title = "Cohort B: genes detected per cell, by specimen", x = NULL, y = "n_feature (log10)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))
ggsave(file.path(fig_dir, "qc_cohort_b_nfeature_violin.png"), p1b, width = 10, height = 5, dpi = 150)

p2 <- ggplot(qc, aes(x = sample_id, y = percent_mt, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE) +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(title = "Cohort B: % mitochondrial counts per cell, by specimen", x = NULL, y = "% MT") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))
ggsave(file.path(fig_dir, "qc_cohort_b_percentmt_violin.png"), p2, width = 10, height = 5, dpi = 150)

p3 <- ggplot(qc, aes(x = n_count, y = n_feature, color = condition)) +
  geom_point(alpha = 0.1, size = 0.3) +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = CONDITION_COLORS) +
  labs(title = "Cohort B: n_count vs n_feature per cell", x = "n_count (log10)", y = "n_feature (log10)") +
  theme_minimal()
ggsave(file.path(fig_dir, "qc_cohort_b_count_vs_feature.png"), p3, width = 7, height = 6, dpi = 150)

p3b <- linear_count_vs_feature_plot(qc, "Cohort B: n_count vs n_feature per cell (linear axes)")
ggsave(file.path(fig_dir, "qc_cohort_b_count_vs_feature_linear.png"), p3b, width = 7, height = 6, dpi = 150)

## B-specific extras -----------------------------------------------------------
p4 <- ggplot(qc, aes(x = TotalUMI, y = n_count)) +
  geom_point(alpha = 0.1, size = 0.3) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Cohort B: recomputed n_count vs supplied TotalUMI",
       x = "supplied TotalUMI (log10)", y = "recomputed n_count (log10)") +
  theme_minimal()
ggsave(file.path(fig_dir, "qc_cohort_b_count_consistency.png"), p4, width = 6, height = 6, dpi = 150)

p5 <- ggplot(qc, aes(x = MT.ratio_pct, y = percent_mt)) +
  geom_point(alpha = 0.1, size = 0.3) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Cohort B: recomputed %MT vs supplied MT.ratio",
       x = "supplied MT.ratio (%)", y = "recomputed percent_mt (%)") +
  theme_minimal()
ggsave(file.path(fig_dir, "qc_cohort_b_mt_consistency.png"), p5, width = 6, height = 6, dpi = 150)

p6 <- ggplot(qc, aes(x = sample_id, y = doublest.score, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE) +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(title = "Cohort B: doublet score per specimen (supplied)", x = NULL, y = "doublest.score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))
ggsave(file.path(fig_dir, "qc_cohort_b_doublet_violin.png"), p6, width = 10, height = 5, dpi = 150)

p7 <- ggplot(per_specimen, aes(x = reorder(sample_id, n_cells), y = n_cells, fill = condition)) +
  geom_col() +
  scale_fill_manual(values = CONDITION_COLORS) +
  coord_flip() +
  labs(title = "Cohort B: cells per specimen", x = NULL, y = "n_cells") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6))
ggsave(file.path(fig_dir, "qc_cohort_b_cells_per_specimen.png"), p7, width = 7, height = 8, dpi = 150)

cat("\nWrote 8 figures (5 standard + 3 B-specific) and 2 tables.\n")
