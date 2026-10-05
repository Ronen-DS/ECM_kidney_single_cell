## Sensitivity check (read-only, metadata-only): same retention logic as
## qc_population_retention.R, with cohort A's percent_mt_max raised from the
## applied 2% to 5%, to quantify how much of the disproportionate immune-
## population loss at 2% is specifically an MT-threshold artifact. This does
## NOT change the actual applied QC (still 2% everywhere else in the
## pipeline) - it is a what-if comparison only. Cohort B is included
## unchanged (10%, as applied) purely so the shared build_meta()/aggregate
## logic runs identically; only cohort A's number differs from the original.
##
## Thresholds used here:
##   A: nFeature_RNA > 500 & < 4500, nCount_RNA > 1000, percent_mt < 5  (vs 2 applied)
##   B: nFeature_RNA > 500 & < 4500, nCount_RNA > 1000, percent_mt < 10 (unchanged, as applied)
##
## Output: results/qc_population_retention_A5pct/*.tsv
## Run: Rscript R/qc_population_retention_A5pct.R

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
source(file.path(repo, "R", "integration_utils.R"))  # add_major_population()
library(Seurat)

out_dir <- file.path(repo, "results", "qc_population_retention_A5pct")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

label_mapping <- read.delim(LABEL_MAPPING_PATH, stringsAsFactors = FALSE)

## Only A's percent_mt_max changed (2 -> 5); B kept at its actually-applied 10.
thresholds <- list(
  A = list(n_feature_min = 500, n_feature_max = 4500, n_count_min = 1000, n_count_max = Inf, percent_mt_max = 5),
  B = list(n_feature_min = 500, n_feature_max = 4500, n_count_min = 1000, n_count_max = Inf, percent_mt_max = 10)
)
MIN_SPECIMEN_CELLS <- 20

build_meta <- function(cohort) {
  cat(sprintf("Loading seurat_%s_raw.rds...\n", cohort))
  obj <- readRDS(file.path(SCRATCH_DIR, sprintf("seurat_%s_raw.rds", cohort)))
  obj <- add_major_population(obj, label_mapping)
  meta <- obj@meta.data
  meta$cell_key <- rownames(meta)
  rm(obj); gc(verbose = FALSE)

  th <- thresholds[[cohort]]
  meta$fail_nfeature_low  <- meta$nFeature_RNA <= th$n_feature_min
  meta$fail_nfeature_high <- meta$nFeature_RNA >= th$n_feature_max
  meta$fail_ncount_low    <- meta$nCount_RNA <= th$n_count_min
  meta$fail_percent_mt    <- meta$percent_mt >= th$percent_mt_max
  meta$cell_pass_qc <- !(meta$fail_nfeature_low | meta$fail_nfeature_high |
                          meta$fail_ncount_low | meta$fail_percent_mt)

  post_filter_counts <- table(meta$sample_id[meta$cell_pass_qc])
  keep_specimens <- names(post_filter_counts)[post_filter_counts > MIN_SPECIMEN_CELLS]
  meta$specimen_dropped <- !(meta$sample_id %in% keep_specimens)
  meta$pass_final <- meta$cell_pass_qc & !meta$specimen_dropped
  meta
}

meta_a <- build_meta("A")
meta_b <- build_meta("B")
all_meta <- rbind(
  meta_a[, c("cell_key","sample_id","donor_id","condition","cohort","major_population",
             "fail_nfeature_low","fail_nfeature_high","fail_ncount_low","fail_percent_mt",
             "cell_pass_qc","specimen_dropped","pass_final")],
  meta_b[, c("cell_key","sample_id","donor_id","condition","cohort","major_population",
             "fail_nfeature_low","fail_nfeature_high","fail_ncount_low","fail_percent_mt",
             "cell_pass_qc","specimen_dropped","pass_final")]
)

## --- (a) Retention table: cohort x population x condition ----------------
agg_before <- aggregate(cell_key ~ cohort + major_population + condition, data = all_meta, FUN = length)
names(agg_before)[4] <- "n_before"
agg_after <- aggregate(cell_key ~ cohort + major_population + condition, data = all_meta[all_meta$pass_final, ],
                        FUN = length)
names(agg_after)[4] <- "n_after"

retention <- merge(agg_before, agg_after, by = c("cohort","major_population","condition"), all.x = TRUE)
retention$n_after[is.na(retention$n_after)] <- 0
retention$fraction_retained <- retention$n_after / retention$n_before

cohort_overall <- aggregate(cbind(n_before = n_before, n_after = n_after) ~ cohort, data = retention, FUN = sum)
cohort_overall$overall_fraction_retained <- cohort_overall$n_after / cohort_overall$n_before
retention <- merge(retention, cohort_overall[, c("cohort","overall_fraction_retained")], by = "cohort")
retention$pp_below_cohort_overall <- 100 * (retention$overall_fraction_retained - retention$fraction_retained)
retention$flagged <- retention$pp_below_cohort_overall > 15

retention <- retention[order(retention$cohort, -retention$pp_below_cohort_overall), ]
write.table(retention, file.path(out_dir, "retention_by_population_condition.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== Overall retention per cohort (A at 5% MT) ===\n")
print(cohort_overall)

cat(sprintf("\n=== Flagged groups (>15pp below their cohort's overall retention) ===\n"))
flagged <- retention[retention$flagged, ]
print(flagged[, c("cohort","major_population","condition","n_before","n_after","fraction_retained","pp_below_cohort_overall")],
      row.names = FALSE)

## --- Endothelial paired-donor coverage, before vs after, at 5% MT --------
focal_pops <- c("Endothelial")
coverage_rows <- list()
for (pop in focal_pops) {
  for (stage in c("before", "after")) {
    sub <- all_meta[all_meta$major_population == pop, ]
    if (stage == "after") sub <- sub[sub$pass_final, ]
    if (nrow(sub) == 0) next
    counts <- aggregate(cell_key ~ donor_id + cohort + condition, data = sub, FUN = length)
    names(counts)[4] <- "n_cells"
    wide <- reshape(counts, idvar = c("donor_id","cohort"), timevar = "condition", direction = "wide")
    names(wide) <- sub("^n_cells\\.", "", names(wide))
    if (!"Tumor" %in% names(wide)) wide$Tumor <- 0
    if (!"Adjacent" %in% names(wide)) wide$Adjacent <- 0
    wide$Tumor[is.na(wide$Tumor)] <- 0
    wide$Adjacent[is.na(wide$Adjacent)] <- 0
    wide$population <- pop
    wide$stage <- stage
    coverage_rows[[paste(pop, stage)]] <- wide[, c("population","stage","donor_id","cohort","Tumor","Adjacent")]
  }
}
coverage_tbl <- do.call(rbind, coverage_rows)
write.table(coverage_tbl, file.path(out_dir, "endothelial_paired_donor_coverage_before_after.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== Endothelial paired-donor coverage (>=20 cells both conditions), before vs after QC (A at 5% MT) ===\n")
for (stage in c("before", "after")) {
  sub <- coverage_tbl[coverage_tbl$stage == stage, ]
  qualifying <- sub[sub$Tumor >= 20 & sub$Adjacent >= 20, ]
  cat(sprintf("\nEndothelial (%s), >= 20 both conditions: %d donor(s)\n", stage, nrow(qualifying)))
  if (nrow(qualifying) > 0) print(qualifying[, c("donor_id","cohort","Tumor","Adjacent")], row.names = FALSE)
}

## --- One-line summary numbers for the assignment report -------------------
a_overall <- cohort_overall$overall_fraction_retained[cohort_overall$cohort == "A"]
a_lymphoid_adj <- retention[retention$cohort == "A" & retention$major_population == "Lymphoid" & retention$condition == "Adjacent", ]
cat(sprintf("\n=== Report line numbers ===\nA overall retention at 5%%: %.1f%%\n", 100 * a_overall))
if (nrow(a_lymphoid_adj) == 1) {
  cat(sprintf("A Adjacent Lymphoid retention at 5%%: %.1f%% (%d -> %d cells)\n",
              100 * a_lymphoid_adj$fraction_retained, a_lymphoid_adj$n_before, a_lymphoid_adj$n_after))
} else {
  cat("A Adjacent Lymphoid row not found or not unique - check retention_by_population_condition.tsv directly.\n")
}

cat("\nDone. Tables written to results/qc_population_retention_A5pct/\n")
