## Cohort C QC: recomputes nCount/nFeature/%MT per cell. No author QC fields
## or cell labels exist for this cohort - these diagnostics are the ONLY
## available signal, so treated as more provisional/exploratory here than
## for A or B (nothing to validate them against). Fast enough (H5 reads) to
## just load all 13 specimens directly rather than stream-and-discard.
##
## Run: Rscript R/qc_cohort_c.R

.args <- commandArgs(trailingOnly = FALSE)
.script_path <- sub("--file=", "", .args[grep("--file=", .args)])
repo <- normalizePath(file.path(dirname(.script_path), ".."))
readRenviron(file.path(repo, ".Renviron"))
source(file.path(repo, "config.R"))
source(file.path(repo, "R", "load_cohort_c.R"))
source(file.path(repo, "R", "qc_utils.R"))
library(ggplot2)

sample_manifest <- read.delim(SAMPLE_MANIFEST_PATH, stringsAsFactors = FALSE)
c_manifest <- sample_manifest[sample_manifest$cohort == "C", , drop = FALSE]
c_samples_meta <- read.delim(COHORT_C_SAMPLES_PATH, stringsAsFactors = FALSE)

qc_list <- vector("list", nrow(c_manifest))
for (i in seq_len(nrow(c_manifest))) {
  sample_id <- c_manifest$sample_id[i]
  cat(sprintf("[C QC %d/%d] %s\n", i, nrow(c_manifest), sample_id))
  h5_path <- file.path(DATA_ROOT, c_manifest$file[i])
  res <- load_cohort_c_sample(sample_id, h5_path)

  mt_ids <- mt_row_ids(res$features, id_col = "gene_id")
  m <- compute_qc_metrics(res$matrix, mt_ids)
  m$sample_id <- sample_id
  qc_list[[i]] <- m
  rm(res)
}
qc <- do.call(rbind, qc_list)
qc <- merge(qc, c_manifest[, c("sample_id", "donor_id", "condition")], by = "sample_id")

cat("\nTotal cells:", nrow(qc), "\n")

## --- Per-specimen summary table ---------------------------------------------
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
per_specimen$low_recovery_flag <- per_specimen$n_cells < LOW_RECOVERY_THRESHOLD
write.table(per_specimen, file.path(repo, "results", "tables", "qc_cohort_c_per_specimen.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("\nSpecimens flagged for low cell recovery (<%d):\n", LOW_RECOVERY_THRESHOLD))
print(per_specimen[per_specimen$low_recovery_flag, c("sample_id", "donor_id", "condition", "n_cells")],
      row.names = FALSE)

## --- Figures -----------------------------------------------------------------
fig_dir <- file.path(repo, "results", "figures", "C")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

qc$condition <- factor(qc$condition, levels = c("Tumor", "Adjacent"))
specimen_order <- per_specimen$sample_id[order(per_specimen$condition, per_specimen$n_cells)]
qc$sample_id <- factor(qc$sample_id, levels = specimen_order)
per_specimen$condition <- factor(per_specimen$condition, levels = c("Tumor", "Adjacent"))

p1 <- ggplot(qc, aes(x = sample_id, y = n_count, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE) +
  scale_y_log10() +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(title = "Cohort C: UMI counts per cell, by specimen (author QC criteria not provided)",
       x = NULL, y = "n_count (log10)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7))
ggsave(file.path(fig_dir, "qc_cohort_c_ncount_violin.png"), p1, width = 9, height = 5, dpi = 150)

p2 <- ggplot(qc, aes(x = sample_id, y = percent_mt, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE) +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(title = "Cohort C: % mitochondrial counts per cell, by specimen",
       x = NULL, y = "% MT") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7))
ggsave(file.path(fig_dir, "qc_cohort_c_percentmt_violin.png"), p2, width = 9, height = 5, dpi = 150)

p2b <- ggplot(qc, aes(x = sample_id, y = n_feature, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE) +
  scale_y_log10() +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(title = "Cohort C: genes detected per cell, by specimen", x = NULL, y = "n_feature (log10)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7))
ggsave(file.path(fig_dir, "qc_cohort_c_nfeature_violin.png"), p2b, width = 9, height = 5, dpi = 150)

p3 <- ggplot(qc, aes(x = n_count, y = n_feature, color = condition)) +
  geom_point(alpha = 0.2, size = 0.5) +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = CONDITION_COLORS) +
  labs(title = "Cohort C: n_count vs n_feature per cell", x = "n_count (log10)", y = "n_feature (log10)") +
  theme_minimal()
ggsave(file.path(fig_dir, "qc_cohort_c_count_vs_feature.png"), p3, width = 7, height = 6, dpi = 150)

p3b <- linear_count_vs_feature_plot(qc, "Cohort C: n_count vs n_feature per cell (linear axes)")
ggsave(file.path(fig_dir, "qc_cohort_c_count_vs_feature_linear.png"), p3b, width = 7, height = 6, dpi = 150)

p4 <- ggplot(per_specimen, aes(x = reorder(sample_id, n_cells), y = n_cells, fill = condition)) +
  geom_col() +
  scale_fill_manual(values = CONDITION_COLORS) +
  coord_flip() +
  labs(title = "Cohort C: cells per specimen", x = NULL, y = "n_cells") +
  theme_minimal()
ggsave(file.path(fig_dir, "qc_cohort_c_cells_per_specimen.png"), p4, width = 7, height = 6, dpi = 150)

cat("\nWrote 6 figures and 1 table.\n")
