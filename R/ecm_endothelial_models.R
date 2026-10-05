## Step 2-6 of the Endothelial ECM pseudobulk Tumor-vs-Adjacent comparison.
## Loads the pseudobulk checkpoint from ecm_endothelial_pseudobulk_build.R and
## runs 5 limma-voom models (+ edgeR normalization), panel-level summaries,
## and donor-level paired plots. No single-cell/integrated-assay data is
## touched here - everything is pseudobulk raw counts -> edgeR TMM -> voom.
##
## Gene universe per model: the INTERSECTION of gene symbols measured in the
## cohorts that model actually uses (A alone for model a; A n B for model b;
## A n B n C for c/d/e) - NOT a union with zero-fill. A gene absent from one
## cohort's reference is not a measured zero for that cohort (see the
## PECAM1-in-B precedent established earlier in this project) and must not be
## silently treated as one just because it's being compared across cohorts.
##
## Condition is coded Adjacent=reference, Tumor=treatment throughout, so
## every reported logFC is "Tumor vs Adjacent" (positive = higher in Tumor).
##
## Output (results/ecm_endothelial/):
##   model_{a_A_only,b_AB,c_ABC,d_paired,e_interaction}_panel_genes.tsv
##   model_*_genes_tested.tsv (n genes in/out of filterByExpr, panel gene status)
##   panel_summary.tsv (per model: mean logFC, summed-panel-feature logFC/CI/P)
##   donor_level_diff.tsv (paired donors: per-gene & panel-sum Tumor-Adjacent logCPM diff)
##   plots/donor_diff_strip_<gene>.png (x8: 10 genes would be a lot - one combined faceted plot instead)
##   plots/paired_lines_SPARC_FN1_panelsum.png
##
## Run: Rscript R/ecm_endothelial_models.R

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
user_lib <- file.path(Sys.getenv("USERPROFILE"), "R", "library-x64")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
library(Matrix)
library(edgeR)
library(limma)
library(ggplot2)

out_dir <- file.path(repo, "results", "ecm_endothelial")
plot_dir <- file.path(out_dir, "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

ecm_genes <- c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1",
               "COL6A2", "FN1", "SPARC", "DCN", "LUM")

cat("Loading pseudobulk_endothelial.rds...\n")
pb <- readRDS(file.path(SCRATCH_DIR, "pseudobulk_endothelial.rds"))
pb_counts <- pb$counts
samples_all <- pb$sample_table
samples_all$condition <- factor(samples_all$condition, levels = c("Adjacent", "Tumor"))

## --- Helpers ---------------------------------------------------------------

#' Build the gene-intersected count matrix + sample table for a set of
#' cohorts, optionally restricted to donors with BOTH conditions present
#' (paired_only). Returns y_full (TMM-normalized, ALL intersected genes - used
#' for lib.size/norm.factors and for descriptive/panel-sum work) and
#' y (filterByExpr-retained subset, used for the actual voom fit).
build_dge <- function(cohorts, paired_only = FALSE) {
  samples <- samples_all[samples_all$cohort %in% cohorts, ]
  if (paired_only) {
    key <- paste(samples$cohort, samples$donor_id)
    paired_keys <- names(table(key))[table(key) == 2]
    samples <- samples[key %in% paired_keys, ]
  }
  samples <- samples[order(samples$cohort, samples$donor_id, samples$condition), ]
  genes <- Reduce(intersect, lapply(cohorts, function(co) rownames(pb_counts[[co]])))
  mats <- lapply(cohorts, function(co) {
    cols <- samples$sample_name[samples$cohort == co]
    pb_counts[[co]][genes, cols, drop = FALSE]
  })
  counts <- do.call(cbind, mats)
  counts <- counts[, samples$sample_name, drop = FALSE]  # re-order to match samples
  samples$cohort <- factor(samples$cohort, levels = cohorts)
  samples$donor_id <- factor(samples$donor_id)

  counts_dense <- as.matrix(counts)
  cat(sprintf("  [build_dge %s] genes=%d, samples=%d, dense counts dim=%d x %d, sum=%.0f\n",
              paste(cohorts, collapse = "+"), length(genes), nrow(samples),
              nrow(counts_dense), ncol(counts_dense), sum(counts_dense)))
  stopifnot(nrow(counts_dense) > 0, ncol(counts_dense) > 0)

  y_full <- DGEList(counts = counts_dense, group = samples$condition)
  y_full <- calcNormFactors(y_full, method = "TMM")
  keep <- filterByExpr(y_full, group = samples$condition)
  y <- y_full[keep, , keep.lib.sizes = FALSE]

  list(y_full = y_full, y = y, samples = samples, n_genes_universe = length(genes),
       n_genes_kept = sum(keep), n_genes_total = length(keep))
}

#' n donors, paired vs unpaired, within a sample table.
donor_pairing_summary <- function(samples) {
  key <- paste(samples$cohort, samples$donor_id)
  tab <- table(key)
  sprintf("n_donors=%d (paired=%d, unpaired=%d)", length(tab), sum(tab == 2), sum(tab == 1))
}

#' Extracts panel-gene rows from a limma topTable-derived full-gene table,
#' adding rows (with NA stats) for any panel gene missing from the model's
#' gene universe or dropped by filterByExpr (kept for descriptive output only,
#' per task spec).
panel_rows <- function(tt_full, gene_universe, genes_tested, model_label) {
  rows <- list()
  for (g in ecm_genes) {
    if (!(g %in% gene_universe)) {
      rows[[g]] <- data.frame(model = model_label, gene = g, logFC = NA, CI.L = NA, CI.R = NA,
                               P.Value = NA, adj.P.Val = NA, status = "absent from gene universe (not measured in one+ contributing cohort)")
    } else if (!(g %in% genes_tested)) {
      rows[[g]] <- data.frame(model = model_label, gene = g, logFC = NA, CI.L = NA, CI.R = NA,
                               P.Value = NA, adj.P.Val = NA, status = "filtered by filterByExpr (descriptive only)")
    } else {
      r <- tt_full[g, ]
      rows[[g]] <- data.frame(model = model_label, gene = g, logFC = r$logFC, CI.L = r$CI.L, CI.R = r$CI.R,
                               P.Value = r$P.Value, adj.P.Val = r$adj.P.Val, status = "tested")
    }
  }
  do.call(rbind, rows)
}

#' Single-feature test of the summed panel raw counts, offset by the model's
#' own TMM-normalized effective library size (lib.size * norm.factors from
#' y_full, i.e. computed across the WHOLE transcriptome before filterByExpr -
#' standard TMM practice). Reuses the same design/block structure as the
#' parent model so it is directly comparable to that model's per-gene logFCs.
panel_sum_test <- function(dge, design, block = NULL, label) {
  genes_present <- intersect(ecm_genes, rownames(dge$y_full))
  panel_sum <- colSums(dge$y_full$counts[genes_present, , drop = FALSE])
  eff_lib <- dge$y_full$samples$lib.size * dge$y_full$samples$norm.factors
  logcpm <- log2((panel_sum + 0.5) / (eff_lib / 1e6))
  m <- matrix(logcpm, nrow = 1, dimnames = list("panel_sum", names(logcpm)))
  if (!is.null(block)) {
    cf <- duplicateCorrelation(m, design, block = block)
    fit <- lmFit(m, design, block = block, correlation = cf$consensus.correlation)
  } else {
    fit <- lmFit(m, design)
  }
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = "conditionTumor", number = 1, confint = 0.95)
  data.frame(model = label, logFC = tt$logFC, CI.L = tt$CI.L, CI.R = tt$CI.R,
             P.Value = tt$P.Value, n_genes_in_sum = length(genes_present),
             genes_in_sum = paste(genes_present, collapse = ","))
}

run_voom_model <- function(dge, design, block_var, label) {
  v <- voom(dge$y, design)
  if (!is.null(block_var)) {
    corfit1 <- duplicateCorrelation(v, design, block = block_var)
    v <- voom(dge$y, design, block = block_var, correlation = corfit1$consensus.correlation)
    corfit2 <- duplicateCorrelation(v, design, block = block_var)
    fit <- lmFit(v, design, block = block_var, correlation = corfit2$consensus.correlation)
    cat(sprintf("  [%s] duplicateCorrelation consensus = %.3f\n", label, corfit2$consensus.correlation))
  } else {
    fit <- lmFit(v, design)
  }
  fit <- eBayes(fit)
  fit
}

all_panel_tables <- list()
all_panel_summary <- list()
all_genes_tested <- list()

## === Model (a): single-cohort only ("A only", "B only", "C only"), =========
## ~ condition, donor block. Each cohort's OWN full gene universe is used
## (not an intersection with the other cohorts) - same logic as model a in
## the original run, just generalized to all three cohorts individually so
## each cohort's standalone signal can be compared before any pooling.
run_single_cohort_model <- function(cohort_letter, label) {
  cat(sprintf("\n=== Model %s: %s only ===\n", label, cohort_letter))
  dge <- build_dge(cohort_letter, paired_only = FALSE)
  cat(sprintf("Genes: %d in %s universe, %d pass filterByExpr. Panel genes retained: %s\n",
              dge$n_genes_universe, cohort_letter, dge$n_genes_kept,
              paste(ecm_genes[ecm_genes %in% rownames(dge$y)], collapse = ",")))
  missing_panel <- setdiff(ecm_genes, rownames(dge$y))
  if (length(missing_panel) > 0) cat(sprintf("  Panel genes filtered out (descriptive only): %s\n", paste(missing_panel, collapse = ",")))
  cat(sprintf("  %s\n", donor_pairing_summary(dge$samples)))

  design <- model.matrix(~ condition, data = dge$samples)
  fit <- run_voom_model(dge, design, dge$samples$donor_id, label)
  tt <- topTable(fit, coef = "conditionTumor", number = Inf, confint = 0.95)
  panel <- panel_rows(tt, rownames(dge$y_full), rownames(dge$y), label)
  all_panel_tables[[label]] <<- panel
  all_genes_tested[[label]] <<- data.frame(model = label, n_genes_universe = dge$n_genes_universe,
                                            n_genes_tested = dge$n_genes_kept, n_donors = length(unique(dge$samples$donor_id)))

  sum_res <- panel_sum_test(dge, design, dge$samples$donor_id, label)
  all_panel_summary[[label]] <<- data.frame(model = label, mean_gene_logFC = mean(panel$logFC, na.rm = TRUE),
                                             n_genes_in_mean = sum(!is.na(panel$logFC)),
                                             panel_sum_logFC = sum_res$logFC, panel_sum_CI.L = sum_res$CI.L,
                                             panel_sum_CI.R = sum_res$CI.R, panel_sum_P = sum_res$P.Value)
  dge
}

dge_a <- run_single_cohort_model("A", "a_A_only")
dge_b_only <- run_single_cohort_model("B", "a2_B_only")
dge_c_only <- run_single_cohort_model("C", "a3_C_only")

## === Model (b): A+B, ~ cohort + condition, donor block =====================
cat("\n=== Model b: A+B ===\n")
dge_b <- build_dge(c("A", "B"), paired_only = FALSE)
cat(sprintf("Genes: %d in A n B universe, %d pass filterByExpr.\n", dge_b$n_genes_universe, dge_b$n_genes_kept))
missing_panel_b <- setdiff(ecm_genes, rownames(dge_b$y))
if (length(missing_panel_b) > 0) cat(sprintf("  Panel genes filtered out (descriptive only): %s\n", paste(missing_panel_b, collapse = ",")))
cat(sprintf("  %s\n", donor_pairing_summary(dge_b$samples)))

design_b <- model.matrix(~ cohort + condition, data = dge_b$samples)
fit_b <- run_voom_model(dge_b, design_b, dge_b$samples$donor_id, "b_AB")
tt_b <- topTable(fit_b, coef = "conditionTumor", number = Inf, confint = 0.95)
panel_b <- panel_rows(tt_b, rownames(dge_b$y_full), rownames(dge_b$y), "b_AB")
all_panel_tables[["b_AB"]] <- panel_b
all_genes_tested[["b_AB"]] <- data.frame(model = "b_AB", n_genes_universe = dge_b$n_genes_universe,
                                          n_genes_tested = dge_b$n_genes_kept, n_donors = length(unique(paste(dge_b$samples$cohort, dge_b$samples$donor_id))))

sum_b <- panel_sum_test(dge_b, design_b, dge_b$samples$donor_id, "b_AB")
all_panel_summary[["b_AB"]] <- data.frame(model = "b_AB", mean_gene_logFC = mean(panel_b$logFC, na.rm = TRUE),
                                           n_genes_in_mean = sum(!is.na(panel_b$logFC)),
                                           panel_sum_logFC = sum_b$logFC, panel_sum_CI.L = sum_b$CI.L,
                                           panel_sum_CI.R = sum_b$CI.R, panel_sum_P = sum_b$P.Value)

## === Model (c): A+B+C, ~ cohort + condition, donor block ===================
cat("\n=== Model c: A+B+C ===\n")
dge_c <- build_dge(c("A", "B", "C"), paired_only = FALSE)
cat(sprintf("Genes: %d in A n B n C universe, %d pass filterByExpr.\n", dge_c$n_genes_universe, dge_c$n_genes_kept))
missing_panel_c <- setdiff(ecm_genes, rownames(dge_c$y))
if (length(missing_panel_c) > 0) cat(sprintf("  Panel genes filtered out (descriptive only): %s\n", paste(missing_panel_c, collapse = ",")))
cat(sprintf("  %s\n", donor_pairing_summary(dge_c$samples)))

design_c <- model.matrix(~ cohort + condition, data = dge_c$samples)
fit_c <- run_voom_model(dge_c, design_c, dge_c$samples$donor_id, "c_ABC")
tt_c <- topTable(fit_c, coef = "conditionTumor", number = Inf, confint = 0.95)
panel_c <- panel_rows(tt_c, rownames(dge_c$y_full), rownames(dge_c$y), "c_ABC")
all_panel_tables[["c_ABC"]] <- panel_c
all_genes_tested[["c_ABC"]] <- data.frame(model = "c_ABC", n_genes_universe = dge_c$n_genes_universe,
                                           n_genes_tested = dge_c$n_genes_kept, n_donors = length(unique(paste(dge_c$samples$cohort, dge_c$samples$donor_id))))

sum_c <- panel_sum_test(dge_c, design_c, dge_c$samples$donor_id, "c_ABC")
all_panel_summary[["c_ABC"]] <- data.frame(model = "c_ABC", mean_gene_logFC = mean(panel_c$logFC, na.rm = TRUE),
                                            n_genes_in_mean = sum(!is.na(panel_c$logFC)),
                                            panel_sum_logFC = sum_c$logFC, panel_sum_CI.L = sum_c$CI.L,
                                            panel_sum_CI.R = sum_c$CI.R, panel_sum_P = sum_c$P.Value)

## === Model (d): paired donors only, A+B+C, ~ donor + condition =============
cat("\n=== Model d: paired donors only, A+B+C ===\n")
dge_d <- build_dge(c("A", "B", "C"), paired_only = TRUE)
cat(sprintf("Genes: %d in A n B n C universe, %d pass filterByExpr.\n", dge_d$n_genes_universe, dge_d$n_genes_kept))
missing_panel_d <- setdiff(ecm_genes, rownames(dge_d$y))
if (length(missing_panel_d) > 0) cat(sprintf("  Panel genes filtered out (descriptive only): %s\n", paste(missing_panel_d, collapse = ",")))
cat(sprintf("  %s\n", donor_pairing_summary(dge_d$samples)))

design_d <- model.matrix(~ donor_id + condition, data = dge_d$samples)
v_d <- voom(dge_d$y, design_d)
fit_d <- lmFit(v_d, design_d)
fit_d <- eBayes(fit_d)
tt_d <- topTable(fit_d, coef = "conditionTumor", number = Inf, confint = 0.95)
panel_d <- panel_rows(tt_d, rownames(dge_d$y_full), rownames(dge_d$y), "d_paired")
all_panel_tables[["d_paired"]] <- panel_d
all_genes_tested[["d_paired"]] <- data.frame(model = "d_paired", n_genes_universe = dge_d$n_genes_universe,
                                              n_genes_tested = dge_d$n_genes_kept, n_donors = length(unique(paste(dge_d$samples$cohort, dge_d$samples$donor_id))))

sum_d <- panel_sum_test(dge_d, design_d, block = NULL, "d_paired")
all_panel_summary[["d_paired"]] <- data.frame(model = "d_paired", mean_gene_logFC = mean(panel_d$logFC, na.rm = TRUE),
                                               n_genes_in_mean = sum(!is.na(panel_d$logFC)),
                                               panel_sum_logFC = sum_d$logFC, panel_sum_CI.L = sum_d$CI.L,
                                               panel_sum_CI.R = sum_d$CI.R, panel_sum_P = sum_d$P.Value)

## === Model (e): A+B+C, cohort:condition interaction, F-test ================
cat("\n=== Model e: A+B+C, cohort x condition interaction ===\n")
design_e <- model.matrix(~ cohort * condition, data = dge_c$samples)
fit_e <- run_voom_model(dge_c, design_e, dge_c$samples$donor_id, "e_interaction")
interaction_coefs <- grep(":", colnames(design_e), value = TRUE)
cat(sprintf("  Interaction coefficients tested jointly (F-test): %s\n", paste(interaction_coefs, collapse = ", ")))
tt_e <- topTable(fit_e, coef = interaction_coefs, number = Inf)  # F-test across coefs, BH-adjusted across all genes
panel_e <- data.frame(model = "e_interaction", gene = ecm_genes)
panel_e$F <- NA; panel_e$P.Value <- NA; panel_e$adj.P.Val <- NA; panel_e$status <- NA
for (i in seq_along(ecm_genes)) {
  g <- ecm_genes[i]
  if (!(g %in% rownames(dge_c$y_full))) { panel_e$status[i] <- "absent from gene universe"; next }
  if (!(g %in% rownames(dge_c$y))) { panel_e$status[i] <- "filtered by filterByExpr (descriptive only)"; next }
  r <- tt_e[g, ]
  panel_e$F[i] <- r$F; panel_e$P.Value[i] <- r$P.Value; panel_e$adj.P.Val[i] <- r$adj.P.Val
  panel_e$status[i] <- "tested"
}
all_panel_tables[["e_interaction"]] <- panel_e

cat("\n=== Interaction F-test, panel genes ===\n")
print(panel_e, row.names = FALSE)

## --- Save per-model tables --------------------------------------------------
for (nm in names(all_panel_tables)[names(all_panel_tables) != "e_interaction"]) {
  write.table(all_panel_tables[[nm]], file.path(out_dir, sprintf("model_%s_panel_genes.tsv", nm)),
              sep = "\t", row.names = FALSE, quote = FALSE)
}
write.table(panel_e, file.path(out_dir, "model_e_interaction_panel_genes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(do.call(rbind, all_genes_tested), file.path(out_dir, "models_genes_tested_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(do.call(rbind, all_panel_summary), file.path(out_dir, "panel_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== Per-model panel gene logFC (Tumor vs Adjacent) ===\n")
for (nm in c("a_A_only", "a2_B_only", "a3_C_only", "b_AB", "c_ABC", "d_paired")) {
  cat(sprintf("\n-- %s --\n", nm))
  print(all_panel_tables[[nm]][, c("gene", "logFC", "CI.L", "CI.R", "P.Value", "adj.P.Val", "status")], row.names = FALSE)
}
cat("\n=== Panel-level summary ===\n")
print(do.call(rbind, all_panel_summary), row.names = FALSE)

## === Step 5: donor-level variability (paired donors, A+B+C) ================
cat("\n=== Donor-level paired Tumor-Adjacent logCPM differences ===\n")
genes_present_d <- intersect(ecm_genes, rownames(dge_d$y_full))
logcpm_d <- cpm(dge_d$y_full, log = TRUE, prior.count = 0.5)[genes_present_d, , drop = FALSE]
panel_sum_counts_d <- colSums(dge_d$y_full$counts[genes_present_d, , drop = FALSE])
eff_lib_d <- dge_d$y_full$samples$lib.size * dge_d$y_full$samples$norm.factors
panel_sum_logcpm_d <- log2((panel_sum_counts_d + 0.5) / (eff_lib_d / 1e6))

long_rows <- list()
donor_keys <- unique(paste(dge_d$samples$cohort, dge_d$samples$donor_id))
for (dk in donor_keys) {
  idx <- which(paste(dge_d$samples$cohort, dge_d$samples$donor_id) == dk)
  s <- dge_d$samples[idx, ]
  tumor_idx <- idx[s$condition == "Tumor"]; adj_idx <- idx[s$condition == "Adjacent"]
  if (length(tumor_idx) != 1 || length(adj_idx) != 1) next  # should not happen - build_dge(paired_only=TRUE) guarantees this
  co <- as.character(dge_d$samples$cohort[tumor_idx]); donor <- as.character(dge_d$samples$donor_id[tumor_idx])
  for (g in genes_present_d) {
    long_rows[[length(long_rows) + 1]] <- data.frame(
      cohort = co, donor_id = donor, gene = g,
      tumor_logcpm = logcpm_d[g, tumor_idx], adjacent_logcpm = logcpm_d[g, adj_idx],
      diff = logcpm_d[g, tumor_idx] - logcpm_d[g, adj_idx]
    )
  }
  long_rows[[length(long_rows) + 1]] <- data.frame(
    cohort = co, donor_id = donor, gene = "panel_sum",
    tumor_logcpm = panel_sum_logcpm_d[tumor_idx], adjacent_logcpm = panel_sum_logcpm_d[adj_idx],
    diff = panel_sum_logcpm_d[tumor_idx] - panel_sum_logcpm_d[adj_idx]
  )
}
donor_diff <- do.call(rbind, long_rows)
write.table(donor_diff, file.path(out_dir, "donor_level_diff.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
print(donor_diff[donor_diff$gene %in% c("SPARC", "FN1", "panel_sum"), ], row.names = FALSE)

## --- Plot 1: faceted strip plot, one point per paired donor, per gene ------
strip_df <- donor_diff
strip_df$gene <- factor(strip_df$gene, levels = c(genes_present_d, "panel_sum"))
p_strip <- ggplot(strip_df, aes(x = cohort, y = diff, color = cohort)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_jitter(width = 0.15, height = 0, size = 2.2) +
  facet_wrap(~gene, nrow = 2) +
  labs(title = "Endothelial, paired donors: Tumor - Adjacent logCPM difference",
       subtitle = sprintf("n=%d paired donors (A=%d, B=%d, C=%d)", length(donor_keys),
                           sum(grepl("^A ", donor_keys)), sum(grepl("^B ", donor_keys)), sum(grepl("^C ", donor_keys))),
       x = NULL, y = "Tumor - Adjacent (logCPM)") +
  theme_minimal() +
  theme(legend.position = "top")
ggsave(file.path(plot_dir, "donor_diff_strip_all_genes.png"), p_strip, width = 11, height = 6, dpi = 150, bg = "white")

## --- Plot 2: paired line plot for SPARC, FN1, panel sum --------------------
line_rows <- list()
for (g in intersect(c("SPARC", "FN1", "panel_sum"), unique(donor_diff$gene))) {
  sub <- donor_diff[donor_diff$gene == g, ]
  for (i in seq_len(nrow(sub))) {
    line_rows[[length(line_rows) + 1]] <- data.frame(cohort = sub$cohort[i], donor_id = sub$donor_id[i], gene = g,
                                                       condition = "Adjacent", logcpm = sub$adjacent_logcpm[i])
    line_rows[[length(line_rows) + 1]] <- data.frame(cohort = sub$cohort[i], donor_id = sub$donor_id[i], gene = g,
                                                       condition = "Tumor", logcpm = sub$tumor_logcpm[i])
  }
}
line_df <- do.call(rbind, line_rows)
line_df$condition <- factor(line_df$condition, levels = c("Adjacent", "Tumor"))
line_df$gene <- factor(line_df$gene, levels = c("SPARC", "FN1", "panel_sum"))
line_df$donor_key <- paste(line_df$cohort, line_df$donor_id)

p_line <- ggplot(line_df, aes(x = condition, y = logcpm, group = donor_key, color = cohort)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 2) +
  facet_wrap(~gene, nrow = 1, scales = "free_y") +
  labs(title = "Endothelial, paired donors: Adjacent -> Tumor (logCPM)",
       subtitle = "SPARC, FN1, and the 10-gene ECM panel sum",
       x = NULL, y = "logCPM") +
  theme_minimal() +
  theme(legend.position = "top")
ggsave(file.path(plot_dir, "paired_lines_SPARC_FN1_panelsum.png"), p_line, width = 10, height = 5, dpi = 150, bg = "white")

cat("\nWrote tables to results/ecm_endothelial/ and plots to results/ecm_endothelial/plots/\n")
cat("Done.\n")
