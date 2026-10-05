## Adds an "A_only" model to the clean (marker-unflagged) Endothelial
## pseudobulk results, mirroring the already-existing B_only/C_only models in
## ecm_endothelial_clean_models.R exactly (same normalization, same
## ~condition design with donor as a duplicateCorrelation block, same >=20-
## cell sample threshold - cohort A's own full gene universe, no intersection
## with other cohorts). The "before" (all-Endothelial) A_only numbers are NOT
## recomputed - they already exist from the original ecm_endothelial_models.R
## run (panel_summary.tsv model "a_A_only", model_a_A_only_panel_genes.tsv)
## and are pulled directly from those files.
##
## Appends rows only - does not touch any existing row in either output file
## (B_only, C_only, AB, ABC, d_paired untouched; this script never recomputes
## them).
##
## Output (appended):
##   results/ecm_endothelial/clean/clean_models_panel_sum_SPARC_FN1.tsv
##   results/ecm_endothelial/clean/before_after_comparison.tsv
##
## Run: Rscript R/ecm_endothelial_add_Aonly_clean.R

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

ecm_genes <- c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1",
               "COL6A2", "FN1", "SPARC", "DCN", "LUM")

cat("Loading pseudobulk_endothelial_clean.rds...\n")
pb <- readRDS(file.path(SCRATCH_DIR, "pseudobulk_endothelial_clean.rds"))
pb_counts <- pb$counts
samples_all <- pb$sample_table
samples_all$condition <- factor(samples_all$condition, levels = c("Adjacent", "Tumor"))

## --- Same helpers as ecm_endothelial_clean_models.R, unchanged ------------
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
  corfit1 <- duplicateCorrelation(v, design, block = block_var)
  v <- voom(dge$y, design, block = block_var, correlation = corfit1$consensus.correlation)
  corfit2 <- duplicateCorrelation(v, design, block = block_var)
  fit <- lmFit(v, design, block = block_var, correlation = corfit2$consensus.correlation)
  cat(sprintf("  [%s] duplicateCorrelation consensus = %.3f\n", label, corfit2$consensus.correlation))
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

## --- Clean A_only model, same recipe as B_only/C_only ----------------------
cat("\n=== Clean model A_only ===\n")
dge <- build_dge("A", paired_only = FALSE)
cat(sprintf("Genes: %d universe, %d pass filterByExpr. %s\n",
            dge$n_genes_universe, dge$n_genes_kept, donor_pairing_summary(dge$samples)))

design <- model.matrix(~ condition, data = dge$samples)
fit <- run_voom_model(dge, design, dge$samples$donor_id, "A_only")
tt <- topTable(fit, coef = "conditionTumor", number = Inf, confint = 0.95)

panel_sum_res <- panel_sum_test(dge, design, dge$samples$donor_id)
sparc_res <- extract_gene(tt, rownames(dge$y_full), rownames(dge$y), "SPARC")
fn1_res <- extract_gene(tt, rownames(dge$y_full), rownames(dge$y), "FN1")

n_donors <- length(unique(paste(dge$samples$cohort, dge$samples$donor_id)))
pairing <- donor_pairing_summary(dge$samples)

new_clean_rows <- rbind(
  data.frame(model = "A_only", feature = "panel_sum", n_donors = n_donors, pairing = pairing,
             logFC = panel_sum_res$logFC, CI.L = panel_sum_res$CI.L, CI.R = panel_sum_res$CI.R,
             P.Value = panel_sum_res$P.Value, status = "tested"),
  data.frame(model = "A_only", feature = "SPARC", n_donors = n_donors, pairing = pairing,
             logFC = sparc_res$logFC, CI.L = sparc_res$CI.L, CI.R = sparc_res$CI.R,
             P.Value = sparc_res$P.Value, status = sparc_res$status),
  data.frame(model = "A_only", feature = "FN1", n_donors = n_donors, pairing = pairing,
             logFC = fn1_res$logFC, CI.L = fn1_res$CI.L, CI.R = fn1_res$CI.R,
             P.Value = fn1_res$P.Value, status = fn1_res$status)
)

cat("\n=== New clean A_only rows ===\n")
print(new_clean_rows, row.names = FALSE)

## --- Append to clean_models_panel_sum_SPARC_FN1.tsv (existing rows untouched) ---
clean_path <- file.path(clean_dir, "clean_models_panel_sum_SPARC_FN1.tsv")
existing_clean <- read.delim(clean_path, stringsAsFactors = FALSE)
stopifnot(!("A_only" %in% existing_clean$model))  # guard against double-append on rerun
updated_clean <- rbind(existing_clean, new_clean_rows)
write.table(updated_clean, clean_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("\nAppended 3 rows to %s (now %d rows)\n", clean_path, nrow(updated_clean)))

## --- Before (all-Endothelial) A_only, pulled directly from existing files ---
orig_panel_summary <- read.delim(file.path(out_dir, "panel_summary.tsv"), stringsAsFactors = FALSE)
ps <- orig_panel_summary[orig_panel_summary$model == "a_A_only", ]
before_panel_sum <- data.frame(model = "A_only", feature = "panel_sum", logFC = ps$panel_sum_logFC,
                                CI.L = ps$panel_sum_CI.L, CI.R = ps$panel_sum_CI.R, P.Value = ps$panel_sum_P,
                                state = "all_Endothelial")

orig_gene_tbl <- read.delim(file.path(out_dir, "model_a_A_only_panel_genes.tsv"), stringsAsFactors = FALSE)
before_gene_rows <- list()
for (g in c("SPARC", "FN1")) {
  gr <- orig_gene_tbl[orig_gene_tbl$gene == g, ]
  before_gene_rows[[g]] <- data.frame(model = "A_only", feature = g, logFC = gr$logFC, CI.L = gr$CI.L,
                                       CI.R = gr$CI.R, P.Value = gr$P.Value, state = "all_Endothelial")
}
before_rows <- rbind(before_panel_sum, do.call(rbind, before_gene_rows))

after_rows <- new_clean_rows[, c("model", "feature", "logFC", "CI.L", "CI.R", "P.Value")]
after_rows$state <- "clean_Endothelial"

new_comparison_rows <- rbind(before_rows, after_rows)
cat("\n=== New before/after A_only rows ===\n")
print(new_comparison_rows, row.names = FALSE)

## --- Append to before_after_comparison.tsv (existing rows untouched) -------
comparison_path <- file.path(clean_dir, "before_after_comparison.tsv")
existing_comparison <- read.delim(comparison_path, stringsAsFactors = FALSE)
stopifnot(!("A_only" %in% existing_comparison$model))
updated_comparison <- rbind(existing_comparison, new_comparison_rows)
updated_comparison <- updated_comparison[order(updated_comparison$feature, updated_comparison$model, updated_comparison$state), ]
write.table(updated_comparison, comparison_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("\nAppended 6 rows to %s (now %d rows)\n", comparison_path, nrow(updated_comparison)))

cat("\nDone.\n")
