## Re-fits the Endothelial ECM pseudobulk models on the CLEAN (mesenchymal-
## marker-unflagged) pseudobulk from ecm_endothelial_mesenchymal_flag.R, using
## the identical normalization/design code as ecm_endothelial_models.R, then
## builds a before/after comparison table (all-Endothelial vs clean-
## Endothelial) for the panel sum, SPARC, and FN1 across models.
##
## Run: Rscript R/ecm_endothelial_clean_models.R

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

out_dir <- file.path(repo, "results", "ecm_endothelial")
clean_dir <- file.path(out_dir, "clean")
dir.create(clean_dir, showWarnings = FALSE, recursive = TRUE)

ecm_genes <- c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1",
               "COL6A2", "FN1", "SPARC", "DCN", "LUM")

cat("Loading pseudobulk_endothelial_clean.rds...\n")
pb <- readRDS(file.path(SCRATCH_DIR, "pseudobulk_endothelial_clean.rds"))
pb_counts <- pb$counts
samples_all <- pb$sample_table
samples_all$condition <- factor(samples_all$condition, levels = c("Adjacent", "Tumor"))

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
  counts <- counts[, samples$sample_name, drop = FALSE]
  samples$cohort <- factor(samples$cohort, levels = cohorts)
  samples$donor_id <- factor(samples$donor_id)

  y_full <- DGEList(counts = as.matrix(counts), group = samples$condition)
  y_full <- calcNormFactors(y_full, method = "TMM")
  keep <- filterByExpr(y_full, group = samples$condition)
  y <- y_full[keep, , keep.lib.sizes = FALSE]
  list(y_full = y_full, y = y, samples = samples, n_genes_universe = length(genes),
       n_genes_kept = sum(keep))
}

donor_pairing_summary <- function(samples) {
  key <- paste(samples$cohort, samples$donor_id)
  tab <- table(key)
  sprintf("n_donors=%d (paired=%d, unpaired=%d)", length(tab), sum(tab == 2), sum(tab == 1))
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
  eBayes(fit)
}

extract_gene <- function(tt_full, gene_universe, genes_tested, gene) {
  if (!(gene %in% gene_universe)) return(data.frame(logFC = NA, CI.L = NA, CI.R = NA, P.Value = NA, status = "absent from gene universe"))
  if (!(gene %in% genes_tested)) return(data.frame(logFC = NA, CI.L = NA, CI.R = NA, P.Value = NA, status = "filtered by filterByExpr"))
  r <- tt_full[gene, ]
  data.frame(logFC = r$logFC, CI.L = r$CI.L, CI.R = r$CI.R, P.Value = r$P.Value, status = "tested")
}

panel_sum_test <- function(dge, design, block = NULL) {
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
  data.frame(logFC = tt$logFC, CI.L = tt$CI.L, CI.R = tt$CI.R, P.Value = tt$P.Value)
}

run_model <- function(label, cohorts, paired_only, fixed_donor = FALSE) {
  cat(sprintf("\n=== Clean model %s ===\n", label))
  dge <- build_dge(cohorts, paired_only = paired_only)
  cat(sprintf("Genes: %d universe, %d pass filterByExpr. %s\n",
              dge$n_genes_universe, dge$n_genes_kept, donor_pairing_summary(dge$samples)))

  if (fixed_donor) {
    design <- model.matrix(~ donor_id + condition, data = dge$samples)
    v <- voom(dge$y, design)
    fit <- eBayes(lmFit(v, design))
    block <- NULL
  } else {
    design <- if (length(cohorts) > 1) model.matrix(~ cohort + condition, data = dge$samples) else model.matrix(~ condition, data = dge$samples)
    fit <- run_voom_model(dge, design, dge$samples$donor_id, label)
    block <- dge$samples$donor_id
  }
  tt <- topTable(fit, coef = "conditionTumor", number = Inf, confint = 0.95)

  panel_sum_res <- panel_sum_test(dge, design, block)
  sparc_res <- extract_gene(tt, rownames(dge$y_full), rownames(dge$y), "SPARC")
  fn1_res <- extract_gene(tt, rownames(dge$y_full), rownames(dge$y), "FN1")

  list(label = label, n_donors = length(unique(paste(dge$samples$cohort, dge$samples$donor_id))),
       pairing = donor_pairing_summary(dge$samples),
       panel_sum = panel_sum_res, SPARC = sparc_res, FN1 = fn1_res)
}

results_clean <- list(
  B_only   = run_model("B_only", "B", paired_only = FALSE),
  C_only   = run_model("C_only", "C", paired_only = FALSE),
  AB       = run_model("AB", c("A", "B"), paired_only = FALSE),
  ABC      = run_model("ABC", c("A", "B", "C"), paired_only = FALSE),
  d_paired = run_model("d_paired", c("A", "B", "C"), paired_only = TRUE, fixed_donor = TRUE)
)

## --- Save clean per-model results -------------------------------------------
rows <- list()
for (nm in names(results_clean)) {
  r <- results_clean[[nm]]
  for (feat in c("panel_sum", "SPARC", "FN1")) {
    v <- r[[feat]]
    rows[[length(rows) + 1]] <- data.frame(model = nm, feature = feat, n_donors = r$n_donors, pairing = r$pairing,
                                            logFC = v$logFC, CI.L = v$CI.L, CI.R = v$CI.R, P.Value = v$P.Value,
                                            status = if (!is.null(v$status)) v$status else "tested")
  }
}
clean_tbl <- do.call(rbind, rows)
write.table(clean_tbl, file.path(clean_dir, "clean_models_panel_sum_SPARC_FN1.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n=== Clean-Endothelial results: panel sum, SPARC, FN1 ===\n")
print(clean_tbl, row.names = FALSE)

## --- Before/after comparison table ------------------------------------------
## "Before" = original all-Endothelial results already saved by
## ecm_endothelial_models.R: panel_summary.tsv (panel_sum) and each
## model_*_panel_genes.tsv (SPARC, FN1 individually).
orig_panel_summary <- read.delim(file.path(out_dir, "panel_summary.tsv"), stringsAsFactors = FALSE)
model_file_map <- c(B_only = "a2_B_only", C_only = "a3_C_only", AB = "b_AB", ABC = "c_ABC", d_paired = "d_paired")

before_rows <- list()
for (nm in names(model_file_map)) {
  orig_label <- model_file_map[[nm]]
  ps <- orig_panel_summary[orig_panel_summary$model == orig_label, ]
  before_rows[[length(before_rows) + 1]] <- data.frame(
    model = nm, feature = "panel_sum", logFC = ps$panel_sum_logFC, CI.L = ps$panel_sum_CI.L,
    CI.R = ps$panel_sum_CI.R, P.Value = ps$panel_sum_P)

  gene_tbl <- read.delim(file.path(out_dir, sprintf("model_%s_panel_genes.tsv", orig_label)), stringsAsFactors = FALSE)
  for (g in c("SPARC", "FN1")) {
    gr <- gene_tbl[gene_tbl$gene == g, ]
    before_rows[[length(before_rows) + 1]] <- data.frame(
      model = nm, feature = g, logFC = gr$logFC, CI.L = gr$CI.L, CI.R = gr$CI.R, P.Value = gr$P.Value)
  }
}
before_tbl <- do.call(rbind, before_rows)
before_tbl$state <- "all_Endothelial"

after_tbl <- clean_tbl[, c("model", "feature", "logFC", "CI.L", "CI.R", "P.Value")]
after_tbl$state <- "clean_Endothelial"

comparison <- rbind(before_tbl, after_tbl)
comparison <- comparison[order(comparison$feature, comparison$model, comparison$state), ]
write.table(comparison, file.path(clean_dir, "before_after_comparison.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== BEFORE (all-Endothelial) vs AFTER (clean-Endothelial): log2FC (95% CI), Tumor vs Adjacent ===\n")
for (feat in c("panel_sum", "SPARC", "FN1")) {
  cat(sprintf("\n-- %s --\n", feat))
  print(comparison[comparison$feature == feat, c("model", "state", "logFC", "CI.L", "CI.R", "P.Value")], row.names = FALSE)
}

cat("\nDone.\n")
