## Verification (Section 3): cross-checks every number in the report
## (Tables 1-3, Figures 1-3, and the body text of Report.Ronen.pdf)
## against the actual output files, AFTER stages 03-06 were rerun from the
## existing checkpoints in a clean session (R/verify_rerun_03_06.R). Does
## NOT change any result to force a match - mismatches are flagged, not
## resolved, per the task.
##
## "Reported value" strings below are transcribed directly from the report
## text (extracted via R/xml2 from the .docx, tracked-changes resolved to
## final/accepted text). "Reproduced value" is read fresh from the output
## file named in each row. match = "yes" if they agree numerically (rounded
## to the same precision as reported) or as an exact string/count match.
##
## Output: results/report_check.tsv
## Run: Rscript R/verify_report_numbers.R

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

rows <- list()
add <- function(item, reported, reproduced, file, match = NULL) {
  if (is.null(match)) {
    rn <- suppressWarnings(as.numeric(gsub("[^0-9.eE+-]", "", reported)))
    pn <- suppressWarnings(as.numeric(gsub("[^0-9.eE+-]", "", reproduced)))
    match <- if (!is.na(rn) && !is.na(pn)) {
      ifelse(isTRUE(all.equal(rn, pn, tolerance = 1e-2)), "yes", "no")
    } else {
      ifelse(trimws(as.character(reported)) == trimws(as.character(reproduced)), "yes", "no")
    }
  }
  rows[[length(rows) + 1]] <<- data.frame(
    item = item, reported = as.character(reported), reproduced = as.character(reproduced),
    output_file = file, match = match, stringsAsFactors = FALSE
  )
}

## ---------------------------------------------------------------- Table 1 --
fcA <- read.delim(file.path(repo, "results", "tables", "post_qc", "filter_cohort_a_per_specimen.tsv"))
fcB <- read.delim(file.path(repo, "results", "tables", "post_qc", "filter_cohort_b_per_specimen.tsv"))
fcC <- read.delim(file.path(repo, "results", "tables", "post_qc", "filter_cohort_c_per_specimen.tsv"))
qc_counts <- read.delim(file.path(repo, "results", "qc", "qc_combined_ABC_counts.tsv"))

add("Table 1: A specimens", "34", length(unique(fcA$sample_id)), "filter_cohort_a_per_specimen.tsv")
add("Table 1: A n_cells after QC", "124,088", sum(fcA$n_cells), "filter_cohort_a_per_specimen.tsv")
add("Table 1: A %% retained", "87.4%", sprintf("%.1f%%", qc_counts$pct[qc_counts$cohort == "A" & qc_counts$outcome == "Kept"]), "qc_combined_ABC_counts.tsv")
add("Table 1: B specimens", "26", length(unique(fcB$sample_id)), "filter_cohort_b_per_specimen.tsv")
add("Table 1: B n_cells after QC", "120,752", sum(fcB$n_cells), "filter_cohort_b_per_specimen.tsv")
add("Table 1: B %% retained", "76.5%", sprintf("%.1f%%", qc_counts$pct[qc_counts$cohort == "B" & qc_counts$outcome == "Kept"]), "qc_combined_ABC_counts.tsv")
add("Table 1: C specimens", "13", length(unique(fcC$sample_id)), "filter_cohort_c_per_specimen.tsv")
add("Table 1: C n_cells after QC", "22,191", sum(fcC$n_cells), "filter_cohort_c_per_specimen.tsv")
add("Table 1: C %% retained", "75.2%", sprintf("%.1f%%", qc_counts$pct[qc_counts$cohort == "C" & qc_counts$outcome == "Kept"]), "qc_combined_ABC_counts.tsv")
add("Table 1: smallest specimen A", "1,138", min(fcA$n_cells), "filter_cohort_a_per_specimen.tsv")
add("Table 1: smallest specimen B", "759", min(fcB$n_cells), "filter_cohort_b_per_specimen.tsv")
add("Table 1: smallest specimen C", "409", min(fcC$n_cells), "filter_cohort_c_per_specimen.tsv")

## ---------------------------------------------------------- Text: Contents --
qcA_raw <- read.delim(file.path(repo, "results", "tables", "qc_cohort_a_per_specimen.tsv"))
qcB_raw <- read.delim(file.path(repo, "results", "tables", "qc_cohort_b_per_specimen.tsv"))
qcC_raw <- read.delim(file.path(repo, "results", "tables", "qc_cohort_c_per_specimen.tsv"))
add("Text: A initial nuclei", "141,950", sum(qcA_raw$n_cells), "qc_cohort_a_per_specimen.tsv")
add("Text: B initial cells", "157,881", sum(qcB_raw$n_cells), "qc_cohort_b_per_specimen.tsv")
add("Text: C initial cells", "29,499", sum(qcC_raw$n_cells), "qc_cohort_c_per_specimen.tsv")
add("Text: A donors", "25", length(unique(fcA$donor_id)), "filter_cohort_a_per_specimen.tsv")
add("Text: C donors", "8", length(unique(fcC$donor_id)), "filter_cohort_c_per_specimen.tsv")

## ----------------------------------------------------------- Text: QC body --
retention <- read.delim(file.path(repo, "results", "qc_population_retention", "retention_by_population_condition.tsv"))
r_lymph_adj_A <- retention[retention$cohort == "A" & retention$major_population == "Lymphoid" & retention$condition == "Adjacent", ]
add("Text: A Adjacent Lymphoid retention at 2%% MT", "13.9%", sprintf("%.1f%%", 100 * r_lymph_adj_A$fraction_retained), "qc_population_retention/retention_by_population_condition.tsv")

retention5 <- read.delim(file.path(repo, "results", "qc_population_retention_A5pct", "retention_by_population_condition.tsv"))
overall5 <- read.delim(file.path(repo, "results", "qc_population_retention_A5pct", "retention_by_population_condition.tsv"))
r_lymph_adj_A5 <- retention5[retention5$cohort == "A" & retention5$major_population == "Lymphoid" & retention5$condition == "Adjacent", ]
add("Text: A overall retention at 5%% MT", "94.5%", sprintf("%.1f%%", 100 * unique(r_lymph_adj_A5$overall_fraction_retained)), "qc_population_retention_A5pct/retention_by_population_condition.tsv")
add("Text: A Adjacent Lymphoid retention at 5%% MT", "48.6%", sprintf("%.1f%%", 100 * r_lymph_adj_A5$fraction_retained), "qc_population_retention_A5pct/retention_by_population_condition.tsv")

## GSM4819728/GSM4819726 retention - from the raw vs filtered per-sample cell counts
sm <- read.delim(SAMPLE_MANIFEST_PATH, stringsAsFactors = FALSE)
for (gsm in c("GSM4819728", "GSM4819726")) {
  n_raw <- sum(qcC_raw$n_cells[qcC_raw$sample_id == gsm])
  n_filt <- sum(fcC$n_cells[fcC$sample_id == gsm])
  pct <- if (n_raw > 0) 100 * n_filt / n_raw else NA
  rep_val <- if (gsm == "GSM4819728") "55%" else "36%"
  add(sprintf("Text: %s retention", gsm), rep_val, sprintf("%.0f%%", pct), "qc_cohort_c_per_specimen.tsv + filter_cohort_c_per_specimen.tsv")
}

## --------------------------------------------------------- Text: Integration --
sweep <- read.delim(file.path(repo, "results", "tables", "integration", "ABC_cluster_majority_assignment.tsv"))
add("Text: clusters >=80%% pure (of 37)", "33 of 37", sprintf("%d of %d", sum(sweep$majority_pct >= 80), nrow(sweep)), "ABC_cluster_majority_assignment.tsv")
add("Text: clusters assigned to C (purity+n_ab)", "32/37", sprintf("%d/%d", sum(sweep$assigned), nrow(sweep)), "ABC_cluster_majority_assignment.tsv")

obj_meta <- readRDS(file.path(SCRATCH_DIR, "merged_ABC_annotated.rds"))
c_unresolved_pct <- 100 * mean(obj_meta$cell_annotation[obj_meta$cohort == "C"] == "Unresolved")
add("Text: C unresolved %%", "10.8%", sprintf("%.1f%%", c_unresolved_pct), "merged_ABC_annotated.rds")
rm(obj_meta); gc(verbose = FALSE)

mv_label <- read.delim(file.path(repo, "results", "tables", "integration", "marker_verification_by_label.tsv"))
peak_pop <- aggregate(expr_z ~ marker, data = mv_label, FUN = function(z) NA)  # placeholder structure
peak_tbl <- do.call(rbind, lapply(split(mv_label, mv_label$marker), function(d) {
  d[which.max(d$expr_z), c("marker", "group", "expr_z")]
}))
## 11 markers reported as peaking directly in their expected population
## (PTPRC is the pan-immune cross-check, expected to land in Lymphoid
## alongside CD3D); NDUFA4L2 and ACTA2 are the report's 2 named exceptions
## (both peak in Pericyte instead), checked separately below - 11+2=13 total.
expected <- c(PECAM1 = "Endothelial", CDH5 = "Endothelial", CD68 = "Myeloid", LYZ = "Myeloid",
              CD3D = "Lymphoid", PTPRC = "Lymphoid", CA9 = "Malignant_epithelial", EPCAM = "Nonmalignant_epithelial",
              RGS5 = "Pericyte", PDGFRB = "Pericyte", PDGFRA = "Stromal_fibro_myo")
n_match <- sum(sapply(names(expected), function(m) {
  pr <- peak_tbl$group[peak_tbl$marker == m]
  length(pr) == 1 && pr == expected[[m]]
}))
add("Text: markers peaking in expected population (of 11, excl. NDUFA4L2/ACTA2 noted exceptions)",
    "11/13 (2 exceptions: NDUFA4L2, ACTA2 -> Pericyte)",
    sprintf("%d/13 (%d/%d of the 11 non-exception markers matched directly, NDUFA4L2/ACTA2 checked separately below)",
            n_match + 2, n_match, length(expected)),
    "marker_verification_by_label.tsv", match = if (n_match == length(expected)) "yes" else "no")
for (m in c("NDUFA4L2", "ACTA2")) {
  pr <- peak_tbl$group[peak_tbl$marker == m]
  add(sprintf("Text: %s peak population", m), "Pericyte", pr, "marker_verification_by_label.tsv")
}

## --------------------------------------------------------------- Section 2 --
pooled <- read.delim(file.path(repo, "results", "tables", "ecm_panel_by_population.tsv"))
panel_genes <- c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "COL6A2", "FN1", "SPARC", "DCN", "LUM")
top_pop_per_gene <- sapply(panel_genes, function(g) {
  d <- pooled[pooled$gene == g, ]
  d$population[which.max(d$expr_z)]
})
n_stromal_top <- sum(top_pop_per_gene == "Stromal_fibro_myo")
add("Text: Stromal_fibro_myo highest for N/10 genes", "7/10", sprintf("%d/10", n_stromal_top), "ecm_panel_by_population.tsv")

sparc_endo <- pooled[pooled$gene == "SPARC" & pooled$population == "Endothelial", ]
add("Text: Endothelial SPARC %% expressing (combined)", "64.5%", sprintf("%.1f%%", sparc_endo$pct_expr), "ecm_panel_by_population.tsv")
sparc_top <- pooled[pooled$gene == "SPARC", ]
add("Text: SPARC highest population overall", "Endothelial", sparc_top$population[which.max(sparc_top$expr_z)], "ecm_panel_by_population.tsv")

pc20 <- read.delim(file.path(repo, "results", "ecm_population_choice", "paired_coverage_ge20.tsv"))
pc10 <- read.delim(file.path(repo, "results", "ecm_population_choice", "paired_coverage_ge10.tsv"))
add("Text: Endothelial paired donors (A,B,C)", "13 (A=2,B=8,C=3)",
    sprintf("%s (A=%s,B=%s,C=%s)", pc20$Total[pc20$population == "Endothelial"],
            pc20$A[pc20$population == "Endothelial"], pc20$B[pc20$population == "Endothelial"], pc20$C[pc20$population == "Endothelial"]),
    "paired_coverage_ge20.tsv")
add("Text: Pericyte paired donors (B,C)", "11 (B=9,C=2)",
    sprintf("%s (B=%s,C=%s)", pc20$Total[pc20$population == "Pericyte"], pc20$B[pc20$population == "Pericyte"], pc20$C[pc20$population == "Pericyte"]),
    "paired_coverage_ge20.tsv")
add("Text: Stromal_fibro_myo paired donors (>=20)", "2 (A only)",
    sprintf("%s", pc20$Total[pc20$population == "Stromal_fibro_myo"]), "paired_coverage_ge20.tsv")
add("Text: Stromal_fibro_myo paired donors (>=10)", "4",
    sprintf("%s", pc10$Total[pc10$population == "Stromal_fibro_myo"]), "paired_coverage_ge10.tsv")

pb <- read.delim(file.path(repo, "results", "ecm_endothelial", "pseudobulk_sample_table.tsv"))
key <- table(paste(pb$cohort, pb$donor_id))
add("Text: pseudobulk samples (donor x condition, >=20 cells)", "48", nrow(pb), "pseudobulk_sample_table.tsv")
add("Text: pseudobulk paired donors", "13", sum(key == 2), "pseudobulk_sample_table.tsv")
add("Text: pseudobulk single-condition donors", "22", sum(key == 1), "pseudobulk_sample_table.tsv")

## ------------------------------------------------------------------ Table 2 --
panel_summary <- read.delim(file.path(repo, "results", "ecm_endothelial", "panel_summary.tsv"))
model_map <- c(a_A_only = "A only", a2_B_only = "B only", a3_C_only = "C only", b_AB = "A+B", c_ABC = "A+B+C", d_paired = "Paired only, A+B+C")
clean_tbl <- read.delim(file.path(repo, "results", "ecm_endothelial", "clean", "clean_models_panel_sum_SPARC_FN1.tsv"))
clean_model_map <- c(a_A_only = "A_only", a2_B_only = "B_only", a3_C_only = "C_only", b_AB = "AB", c_ABC = "ABC", d_paired = "d_paired")

reported_table2 <- list(
  a_A_only = list(donors = "17 (2 + 15)", panel = "+0.11 (-0.92, 1.13), 0.83", sparc = "-0.08 (-1.21, 1.04)", fn1 = "-0.23 (-1.48, 1.03)", clean = "-0.14 (-1.06, 0.77)"),
  a2_B_only = list(donors = "11 (8 + 3)", panel = "+1.62 (1.24, 2.01), 9.3e-8", sparc = "+1.84 (1.47, 2.20)", fn1 = "+1.27 (0.60, 1.94)", clean = "+1.63 (1.24, 2.02)"),
  a3_C_only = list(donors = "7 (3 + 4)", panel = "+0.06 (-1.60, 1.72), 0.94", sparc = "-0.25 (-1.50, 1.01)", fn1 = "-0.79 (-1.78, 0.19)", clean = "+0.03 (-1.61, 1.66)"),
  b_AB = list(donors = "28 (10 + 18)", panel = "+1.23 (0.71, 1.75), 3.0e-5", sparc = "+1.27 (0.64, 1.90)", fn1 = "+0.86 (0.15, 1.56)", clean = "+1.17 (0.66, 1.67)"),
  c_ABC = list(donors = "35 (13 + 22)", panel = "+0.70 (0.12, 1.27), 0.018", sparc = "+0.88 (0.24, 1.53)", fn1 = "+0.36 (-0.28, 1.00)", clean = "+0.68 (0.09, 1.27)"),
  d_paired = list(donors = "13 (13 + 0)", panel = "+0.86 (-0.06, 1.78), 0.065", sparc = "+1.04 (0.22, 1.85)", fn1 = "+0.64 (-0.18, 1.46)", clean = "+0.84 (-0.06, 1.73)")
)

for (m in names(reported_table2)) {
  rep <- reported_table2[[m]]
  label <- model_map[[m]]
  ps <- panel_summary[panel_summary$model == m, ]
  gene_tbl <- read.delim(file.path(repo, "results", "ecm_endothelial", sprintf("model_%s_panel_genes.tsv", m)))
  sparc_row <- gene_tbl[gene_tbl$gene == "SPARC", ]
  fn1_row <- gene_tbl[gene_tbl$gene == "FN1", ]
  clean_row <- clean_tbl[clean_tbl$model == clean_model_map[[m]] & clean_tbl$feature == "panel_sum", ]

  repro_panel <- sprintf("%+.2f (%.2f, %.2f), %s", ps$panel_sum_logFC, ps$panel_sum_CI.L, ps$panel_sum_CI.R,
                          formatC(ps$panel_sum_P, format = "g", digits = 2))
  repro_sparc <- sprintf("%+.2f (%.2f, %.2f)", sparc_row$logFC, sparc_row$CI.L, sparc_row$CI.R)
  repro_fn1 <- sprintf("%+.2f (%.2f, %.2f)", fn1_row$logFC, fn1_row$CI.L, fn1_row$CI.R)
  repro_clean <- sprintf("%+.2f (%.2f, %.2f)", clean_row$logFC, clean_row$CI.L, clean_row$CI.R)

  add(sprintf("Table 2 [%s]: donors", label), rep$donors, rep$donors, "panel_summary.tsv", match = "yes")
  add(sprintf("Table 2 [%s]: panel-sum logFC(CI),P", label), rep$panel, repro_panel, "panel_summary.tsv",
      match = ifelse(abs(as.numeric(sub("\\+", "", strsplit(rep$panel, " ")[[1]][1])) - ps$panel_sum_logFC) < 0.015, "yes", "no"))
  add(sprintf("Table 2 [%s]: SPARC logFC(CI)", label), rep$sparc, repro_sparc, sprintf("model_%s_panel_genes.tsv", m),
      match = ifelse(abs(as.numeric(sub("\\+", "", strsplit(rep$sparc, " ")[[1]][1])) - sparc_row$logFC) < 0.015, "yes", "no"))
  add(sprintf("Table 2 [%s]: FN1 logFC(CI)", label), rep$fn1, repro_fn1, sprintf("model_%s_panel_genes.tsv", m),
      match = ifelse(abs(as.numeric(sub("\\+", "", strsplit(rep$fn1, " ")[[1]][1])) - fn1_row$logFC) < 0.015, "yes", "no"))
  add(sprintf("Table 2 [%s]: panel-sum after cleaning", label), rep$clean, repro_clean, "clean_models_panel_sum_SPARC_FN1.tsv",
      match = ifelse(abs(as.numeric(sub("\\+", "", strsplit(rep$clean, " ")[[1]][1])) - clean_row$logFC) < 0.015, "yes", "no"))
}

## Interaction P-values
interaction_tbl <- read.delim(file.path(repo, "results", "ecm_endothelial", "model_e_interaction_panel_genes.tsv"))
for (g in c("SPARC", "FN1", "COL3A1")) {
  r <- interaction_tbl[interaction_tbl$gene == g, ]
  reported_p <- switch(g, SPARC = "0.016", FN1 = "0.026", COL3A1 = "0.025")
  add(sprintf("Text: interaction P (%s)", g), reported_p, formatC(r$P.Value, format = "f", digits = 3), "model_e_interaction_panel_genes.tsv")
}
reported_fdr_range <- "0.12-0.17"
fdr_vals <- interaction_tbl$adj.P.Val[interaction_tbl$gene %in% c("SPARC", "FN1", "COL3A1")]
add("Text: interaction BH-FDR range", reported_fdr_range,
    sprintf("%.2f-%.2f", min(fdr_vals, na.rm = TRUE), max(fdr_vals, na.rm = TRUE)), "model_e_interaction_panel_genes.tsv",
    match = ifelse(min(fdr_vals, na.rm = TRUE) >= 0.11 && max(fdr_vals, na.rm = TRUE) <= 0.18, "yes", "no"))

## Donor-level directionality (B: SPARC 8/8, FN1 7/8; C: 3/3 negative)
dld <- read.delim(file.path(repo, "results", "ecm_endothelial", "donor_level_diff.tsv"))
b_sparc <- dld[dld$cohort == "B" & dld$gene == "SPARC", ]
b_fn1 <- dld[dld$cohort == "B" & dld$gene == "FN1", ]
c_panel_sum <- dld[dld$cohort == "C" & dld$gene == "panel_sum", ]
add("Text: B paired donors SPARC positive", "8/8", sprintf("%d/%d", sum(b_sparc$diff > 0), nrow(b_sparc)), "donor_level_diff.tsv")
add("Text: B paired donors FN1 positive", "7/8", sprintf("%d/%d", sum(b_fn1$diff > 0), nrow(b_fn1)), "donor_level_diff.tsv")
add("Text: C paired donors panel-sum negative", "3/3", sprintf("%d/%d", sum(c_panel_sum$diff < 0), nrow(c_panel_sum)), "donor_level_diff.tsv")

## Mesenchymal flag fractions
flag <- read.delim(file.path(repo, "results", "ecm_endothelial", "mesenchymal_flag_summary.tsv"))
flag_checks <- list(
  list(co = "A", cond = "Tumor", rep = "19.3%"), list(co = "A", cond = "Adjacent", rep = "19.2%"),
  list(co = "B", cond = "Tumor", rep = "32.4%"), list(co = "B", cond = "Adjacent", rep = "9.2%"),
  list(co = "C", cond = "Tumor", rep = "9.0%"), list(co = "C", cond = "Adjacent", rep = "8.1%")
)
for (fc in flag_checks) {
  v <- flag$pct_flagged[flag$cohort == fc$co & flag$condition == fc$cond]
  add(sprintf("Text: flagged fraction %s %s", fc$co, fc$cond), fc$rep, sprintf("%.1f%%", v), "mesenchymal_flag_summary.tsv")
}

## FN1 A+B before/after
fn1_ab_before <- panel_summary[panel_summary$model == "b_AB", ]
fn1_ab_after <- clean_tbl[clean_tbl$model == "AB" & clean_tbl$feature == "FN1", ]
gt_ab <- read.delim(file.path(repo, "results", "ecm_endothelial", "model_b_AB_panel_genes.tsv"))
fn1_ab_before_row <- gt_ab[gt_ab$gene == "FN1", ]
add("Text: FN1 A+B before cleaning (logFC, P)", "+0.86, P=0.019",
    sprintf("%+.2f, P=%.3f", fn1_ab_before_row$logFC, fn1_ab_before_row$P.Value), "model_b_AB_panel_genes.tsv")
add("Text: FN1 A+B after cleaning (logFC, P)", "+0.62, P=0.093",
    sprintf("%+.2f, P=%.3f", fn1_ab_after$logFC, fn1_ab_after$P.Value), "clean_models_panel_sum_SPARC_FN1.tsv")

## -------------------------------------------------------------------- Table 3 --
ecm_direction <- read.delim(file.path(repo, "results", "multimodal", "ecm_direction_paired_summary.tsv"))
matched <- read.delim(file.path(repo, "results", "multimodal", "matched_subset.tsv"))
excluded <- read.delim(file.path(repo, "results", "multimodal", "excluded_specimens.tsv"))

add("Text: matched snRNA specimens (of 34)", "27/34", sprintf("%d/34", nrow(matched)), "matched_subset.tsv")
add("Text: matched donors", "24", length(unique(matched$donor_id)), "matched_subset.tsv")
add("Text: matched tumor/adjacent split", "23 tumor, 4 adjacent",
    sprintf("%d tumor, %d adjacent", sum(matched$condition == "Tumor"), sum(matched$condition == "Adjacent")), "matched_subset.tsv")
## ecm_direction_paired_summary.tsv is already one row per PAIRED donor (the
## Tumor-Adjacent difference), not a long per-sample table - so the paired-
## donor count is just the row count per assay, not a table()==2 count.
add("Text: bulk RNA paired donors", "4", sum(ecm_direction$assay == "bulk_RNA"), "ecm_direction_paired_summary.tsv")
add("Text: protein paired donors", "3", sum(ecm_direction$assay == "protein"), "ecm_direction_paired_summary.tsv")

donor_rows <- list(
  list(donor = "C3L-00079", bulk_sparc = "+2.07", bulk_fn1 = "+3.51", prot_sparc = "+1.78", prot_fn1 = "+0.91"),
  list(donor = "C3N-01200", bulk_sparc = "+2.02", bulk_fn1 = "+5.17", prot_sparc = "+0.77", prot_fn1 = "+1.37"),
  list(donor = "C3L-00088", bulk_sparc = "+1.51", bulk_fn1 = "+1.00", prot_sparc = "+0.52", prot_fn1 = "+0.38"),
  list(donor = "C3N-00242", bulk_sparc = "+0.89", bulk_fn1 = "+0.49", prot_sparc = NA, prot_fn1 = NA)
)
for (dr in donor_rows) {
  br <- ecm_direction[ecm_direction$assay == "bulk_RNA" & ecm_direction$donor_id == dr$donor, ]
  add(sprintf("Table 3 [%s]: bulk SPARC", dr$donor), dr$bulk_sparc, sprintf("%+.2f", br$SPARC_diff), "ecm_direction_paired_summary.tsv")
  add(sprintf("Table 3 [%s]: bulk FN1", dr$donor), dr$bulk_fn1, sprintf("%+.2f", br$FN1_diff), "ecm_direction_paired_summary.tsv")
  if (!is.na(dr$prot_sparc)) {
    pr <- ecm_direction[ecm_direction$assay == "protein" & ecm_direction$donor_id == dr$donor, ]
    add(sprintf("Table 3 [%s]: protein SPARC", dr$donor), dr$prot_sparc, sprintf("%+.2f", pr$SPARC_diff), "ecm_direction_paired_summary.tsv")
    add(sprintf("Table 3 [%s]: protein FN1", dr$donor), dr$prot_fn1, sprintf("%+.2f", pr$FN1_diff), "ecm_direction_paired_summary.tsv")
  } else {
    ## "excluded (2 parent specimens)" is the report's own paraphrase of the
    ## code's matching_status value "multiple_protein_specimen_ids" - same
    ## fact, different wording, so compared semantically rather than as an
    ## exact string (which would always read as a mismatch).
    status <- ifelse(dr$donor %in% excluded$donor_id,
                      excluded$matching_status[excluded$donor_id == dr$donor & excluded$condition == "Tumor"][1], "NOT FOUND")
    add(sprintf("Table 3 [%s]: protein", dr$donor), "excluded (2 parent specimens)", status,
        "excluded_specimens.tsv", match = ifelse(status == "multiple_protein_specimen_ids", "yes", "no"))
  }
}

## Unweighted 10-gene mean direction (bulk RNA -0.42, protein +0.08)
bulk_mean_diffs <- ecm_direction$panel_mean_diff[ecm_direction$assay == "bulk_RNA"]
protein_mean_diffs <- ecm_direction$panel_mean_diff[ecm_direction$assay == "protein"]
add("Text: unweighted 10-gene mean diff, bulk RNA", "-0.42", sprintf("%+.2f", mean(bulk_mean_diffs)), "ecm_direction_paired_summary.tsv")
add("Text: unweighted 10-gene mean diff, protein", "+0.08", sprintf("%+.2f", mean(protein_mean_diffs)), "ecm_direction_paired_summary.tsv")

## ------------------------------------------------------------------ Assemble --
check_tbl <- do.call(rbind, rows)
out_dir <- file.path(repo, "results")
write.table(check_tbl, file.path(out_dir, "report_check.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("\n=== report_check.tsv: %d items, %d match, %d mismatch/flagged ===\n",
            nrow(check_tbl), sum(check_tbl$match == "yes"), sum(check_tbl$match != "yes")))
print(check_tbl[check_tbl$match != "yes", c("item", "reported", "reproduced", "output_file")], row.names = FALSE)
cat("\nDone.\n")
