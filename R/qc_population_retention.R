## Read-only check: did per-cell QC filtering disproportionately remove any
## author-labeled major population, in cohorts A and B? Cohort C is excluded
## (no author labels to map to a population - see label_mapping.tsv).
##
## This does NOT rerun any filtering - it replicates the exact per-cell and
## per-specimen logic already applied by filter_cohort_a.R / filter_cohort_b.R
## (see qc_utils.R::filter_seurat_cohort) directly against the metadata of the
## UNFILTERED checkpoints, so each cell's pass/fail reason can be recorded
## individually (filter_seurat_cohort() itself only returns the filtered
## object, not a per-cell reason).
##
## Thresholds (already-applied, not re-derived here):
##   A: nFeature_RNA > 500 & < 4500, nCount_RNA > 1000, percent_mt < 2
##   B: nFeature_RNA > 500 & < 4500, nCount_RNA > 1000, percent_mt < 10
##   both: specimen dropped entirely if its POST-cell-filter cell count <= 20
##
## Input:
##   <SCRATCH_DIR>/seurat_A_raw.rds, seurat_B_raw.rds       (unfiltered, pre-QC)
##   <SCRATCH_DIR>/seurat_A_filtered.rds, seurat_B_filtered.rds (post-QC, for cross-check only)
##   metadata/label_mapping.tsv                              (author label -> major_population)
##
## Output (results/qc_population_retention/):
##   retention_by_population_condition.tsv
##   flagged_group_fail_reasons.tsv
##   paired_donor_coverage_before_after.tsv
##
## Run: Rscript R/qc_population_retention.R

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

out_dir <- file.path(repo, "results", "qc_population_retention")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

label_mapping <- read.delim(LABEL_MAPPING_PATH, stringsAsFactors = FALSE)

## QC thresholds exactly as applied in filter_cohort_a.R / filter_cohort_b.R
thresholds <- list(
  A = list(n_feature_min = 500, n_feature_max = 4500, n_count_min = 1000, n_count_max = Inf, percent_mt_max = 2),
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

  ## Specimen-level drop: checked against POST-cell-filter count per specimen,
  ## exactly as filter_seurat_cohort() does (table() on cell_pass_qc==TRUE cells only).
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

## --- Sanity cross-check against the actual seurat_{A,B}_filtered.rds checkpoints ---
for (cohort in c("A", "B")) {
  filt <- readRDS(file.path(SCRATCH_DIR, sprintf("seurat_%s_filtered.rds", cohort)))
  replicated_keep <- sum(all_meta$cohort == cohort & all_meta$pass_final)
  cat(sprintf("Cross-check %s: filtered.rds has %d cells; replicated pass_final = %d (match: %s)\n",
              cohort, ncol(filt), replicated_keep, ncol(filt) == replicated_keep))
  rm(filt); gc(verbose = FALSE)
}

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

cat("\n=== Overall retention per cohort ===\n")
print(cohort_overall)

cat(sprintf("\n=== Flagged groups (>15pp below their cohort's overall retention) ===\n"))
flagged <- retention[retention$flagged, ]
print(flagged[, c("cohort","major_population","condition","n_before","n_after","fraction_retained","pp_below_cohort_overall")],
      row.names = FALSE)

## --- (b) For flagged groups: which single fail criterion removed the most cells ---
fail_reason_rows <- list()
if (nrow(flagged) > 0) {
  for (i in seq_len(nrow(flagged))) {
    g <- flagged[i, ]
    sub <- all_meta[all_meta$cohort == g$cohort & all_meta$major_population == g$major_population &
                     all_meta$condition == g$condition, ]
    failing <- sub[!sub$pass_final, ]
    reason_counts <- c(
      nfeature_low  = sum(failing$fail_nfeature_low),
      nfeature_high = sum(failing$fail_nfeature_high),
      ncount_low    = sum(failing$fail_ncount_low),
      percent_mt    = sum(failing$fail_percent_mt),
      specimen_dropped_only = sum(failing$specimen_dropped & failing$cell_pass_qc)
    )
    top_reason <- names(reason_counts)[which.max(reason_counts)]
    fail_reason_rows[[length(fail_reason_rows) + 1]] <- data.frame(
      cohort = g$cohort, major_population = g$major_population, condition = g$condition,
      n_failing = nrow(failing),
      nfeature_low = reason_counts["nfeature_low"], nfeature_high = reason_counts["nfeature_high"],
      ncount_low = reason_counts["ncount_low"], percent_mt = reason_counts["percent_mt"],
      specimen_dropped_only = reason_counts["specimen_dropped_only"],
      top_reason = top_reason
    )
  }
}
fail_reasons <- if (length(fail_reason_rows) > 0) do.call(rbind, fail_reason_rows) else
  data.frame(cohort=character(), major_population=character(), condition=character())
write.table(fail_reasons, file.path(out_dir, "flagged_group_fail_reasons.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== Flagged-group fail reasons (counts may overlap - a cell can fail >1 criterion) ===\n")
if (nrow(fail_reasons) > 0) print(fail_reasons, row.names = FALSE) else cat("(no flagged groups)\n")

## --- (c) Paired-donor coverage before vs after, for 3 populations ---------
focal_pops <- c("Stromal_fibro_myo", "Endothelial", "Pericyte")
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
write.table(coverage_tbl, file.path(out_dir, "paired_donor_coverage_before_after.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== Paired-donor coverage (>=10 and >=20 cells in BOTH Tumor and Adjacent), before vs after QC ===\n")
for (pop in focal_pops) {
  for (stage in c("before", "after")) {
    sub <- coverage_tbl[coverage_tbl$population == pop & coverage_tbl$stage == stage, ]
    if (nrow(sub) == 0) {
      cat(sprintf("\n%s (%s): no cells at all\n", pop, stage))
      next
    }
    for (thresh in c(10, 20)) {
      qualifying <- sub[sub$Tumor >= thresh & sub$Adjacent >= thresh, ]
      cat(sprintf("\n%s (%s), >= %d both conditions: %d donor(s)\n", pop, stage, thresh, nrow(qualifying)))
      if (nrow(qualifying) > 0) {
        for (j in seq_len(nrow(qualifying))) {
          cat(sprintf("  %s (%s): Tumor=%d / Adjacent=%d\n",
                      qualifying$donor_id[j], qualifying$cohort[j], qualifying$Tumor[j], qualifying$Adjacent[j]))
        }
      }
    }
  }
}

cat("\nDone. Tables written to results/qc_population_retention/\n")
