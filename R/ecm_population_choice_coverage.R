## Read-only coverage check to justify which population (Stromal_fibro_myo,
## Pericyte, Endothelial) best supports a paired donor-level Tumor-vs-Adjacent
## comparison in the integrated A+B+C object. No QC, integration, clustering,
## labels or thresholds are changed here.
##
## Population labels: `cell_annotation` in merged_ABC_annotated.rds (this is
## the actual column name - A/B from label_mapping.tsv via
## add_major_population(), C from resolution-1 cluster-majority transfer;
## see annotate_c_by_cluster_majority.R. The task referred to this column as
## `pop_final`; no such column exists in this object, `cell_annotation` is
## the one actually written by the annotation pipeline).
##
## Donor grouping: sample_manifest.tsv already gives donor_id + condition per
## specimen for ALL THREE cohorts (not just C), and repeated specimens from
## the same donor+condition (e.g. cohort A's C3N-01200 Tumor, split across 3
## archives; cohort C's Patient SS_2014 Adjacent, split into Cortex + Medulla
## per metadata/original/cohort_C_samples.tsv's `tissue` column) already
## share one donor_id/condition pair - so aggregating merged_ABC_annotated.rds
## cells by (cohort, donor_id, condition) automatically pools repeated
## specimens with no extra join needed. cohort_C_samples.tsv is read here only
## to confirm/document the cortex+medulla case; it contributes no columns used
## downstream, since donor_id/condition already live in sample_manifest.tsv
## and therefore in merged_ABC_annotated.rds's own metadata.
##
## Input:  merged_ABC_annotated.rds, metadata/sample_manifest.tsv,
##         metadata/label_mapping.tsv, metadata/original/cohort_C_samples.tsv
## Output: results/ecm_population_choice/{paired_coverage_ge20.tsv,
##         paired_coverage_ge10.tsv, per_donor_detail.tsv,
##         single_condition_donors.tsv}
##
## Run: Rscript R/ecm_population_choice_coverage.R

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

out_dir <- file.path(repo, "results", "ecm_population_choice")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

candidate_pops <- c("Stromal_fibro_myo", "Pericyte", "Endothelial")
cohorts <- c("A", "B", "C")

cat("Loading merged_ABC_annotated.rds for metadata only...\n")
obj <- readRDS(file.path(SCRATCH_DIR, "merged_ABC_annotated.rds"))
meta <- obj@meta.data
rm(obj); gc(verbose = FALSE)
stopifnot(all(c("cohort", "donor_id", "condition", "cell_annotation") %in% colnames(meta)))

sample_manifest <- read.delim(SAMPLE_MANIFEST_PATH, stringsAsFactors = FALSE)
label_mapping <- read.delim(LABEL_MAPPING_PATH, stringsAsFactors = FALSE)
cohort_c_samples <- read.delim(COHORT_C_SAMPLES_PATH, stringsAsFactors = FALSE)

cat("\nDonor+condition combos backed by >1 specimen (sample_manifest.tsv; these already\n")
cat("collapse into one donor_id x condition group when aggregating cell-level metadata):\n")
dup_check <- aggregate(sample_id ~ cohort + donor_id + condition, data = sample_manifest, FUN = length)
print(dup_check[dup_check$sample_id > 1, ], row.names = FALSE)
cat("\ncohort_C_samples.tsv tissue values for Patient SS_2014 (cortex/medulla check):\n")
print(cohort_c_samples[cohort_c_samples$individual == "Patient SS_2014", c("sample","individual","tissue")], row.names = FALSE)

## --- Per-donor detail: full donor x condition grid per cohort, 0-filled ---
build_detail <- function(pop) {
  rows <- list()
  for (co in cohorts) {
    donors <- unique(meta$donor_id[meta$cohort == co])
    grid <- expand.grid(donor_id = donors, condition = c("Tumor", "Adjacent"), stringsAsFactors = FALSE)
    cells <- meta[meta$cohort == co & meta$cell_annotation == pop, ]
    if (nrow(cells) > 0) {
      cnt <- aggregate(dummy ~ donor_id + condition,
                        data = data.frame(dummy = 1, donor_id = cells$donor_id, condition = cells$condition),
                        FUN = length)
      names(cnt)[3] <- "n_cells"
      g <- merge(grid, cnt, by = c("donor_id", "condition"), all.x = TRUE)
    } else {
      g <- grid
      g$n_cells <- NA_real_
    }
    g$n_cells[is.na(g$n_cells)] <- 0
    g$cohort <- co
    rows[[co]] <- g
  }
  full <- do.call(rbind, rows)
  wide <- reshape(full, idvar = c("donor_id", "cohort"), timevar = "condition", direction = "wide")
  names(wide) <- sub("^n_cells\\.", "", names(wide))
  if (!"Tumor" %in% names(wide)) wide$Tumor <- 0
  if (!"Adjacent" %in% names(wide)) wide$Adjacent <- 0
  wide$Tumor[is.na(wide$Tumor)] <- 0
  wide$Adjacent[is.na(wide$Adjacent)] <- 0
  wide$population <- pop
  wide$qualifies_ge10 <- wide$Tumor >= 10 & wide$Adjacent >= 10
  wide$qualifies_ge20 <- wide$Tumor >= 20 & wide$Adjacent >= 20
  wide[order(wide$cohort, -wide$Tumor - wide$Adjacent), c("population", "cohort", "donor_id", "Tumor", "Adjacent", "qualifies_ge10", "qualifies_ge20")]
}

detail_all <- do.call(rbind, lapply(candidate_pops, build_detail))
write.table(detail_all, file.path(out_dir, "per_donor_detail.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## --- Structural-absence detection -----------------------------------------
## A/B: population has no label_mapping row at all for that cohort -> label
## does not exist there (not just zero cells). C has no label_mapping rows by
## construction (cluster-majority transfer, not author labels) - detect
## "no cells assigned" empirically instead (0 total cells of that population
## actually present in C anywhere, not just failing the donor threshold).
label_exists_ab <- function(co, pop) any(label_mapping$cohort == co & label_mapping$population == pop)
c_has_cells <- function(pop) sum(meta$cohort == "C" & meta$cell_annotation == pop) > 0

cell_label <- function(pop, co, n_donors) {
  if (co %in% c("A", "B") && !label_exists_ab(co, pop)) return(sprintf("no label in %s", co))
  if (co == "C" && !c_has_cells(pop)) return("no cells assigned in C")
  as.character(n_donors)
}

## --- Paired coverage tables (rows = population, cols = A, B, C, Total) ----
build_coverage <- function(thresh_col) {
  rows <- list()
  for (pop in candidate_pops) {
    sub <- detail_all[detail_all$population == pop, ]
    counts <- sapply(cohorts, function(co) sum(sub$cohort == co & sub[[thresh_col]]))
    total <- sum(counts)  # donor IDs are cohort-specific strings - no cross-cohort collision, safe to sum
    rows[[pop]] <- data.frame(
      population = pop,
      A = cell_label(pop, "A", counts["A"]),
      B = cell_label(pop, "B", counts["B"]),
      C = cell_label(pop, "C", counts["C"]),
      Total = total,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

coverage_ge20 <- build_coverage("qualifies_ge20")
coverage_ge10 <- build_coverage("qualifies_ge10")
write.table(coverage_ge20, file.path(out_dir, "paired_coverage_ge20.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(coverage_ge10, file.path(out_dir, "paired_coverage_ge10.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== Paired coverage: donors with >= 20 cells in BOTH Tumor and Adjacent ===\n")
print(coverage_ge20, row.names = FALSE)
cat("\n=== Paired coverage: donors with >= 10 cells in BOTH Tumor and Adjacent ===\n")
print(coverage_ge10, row.names = FALSE)

## --- Single-condition-only donors (Tumor-only / Adjacent-only), per population x cohort ---
single_cond_rows <- list()
for (pop in candidate_pops) {
  sub <- detail_all[detail_all$population == pop, ]
  for (co in cohorts) {
    s <- sub[sub$cohort == co, ]
    if (co %in% c("A", "B") && !label_exists_ab(co, pop)) {
      tumor_only <- NA_integer_; adjacent_only <- NA_integer_; note <- sprintf("no label in %s", co)
    } else if (co == "C" && !c_has_cells(pop)) {
      tumor_only <- NA_integer_; adjacent_only <- NA_integer_; note <- "no cells assigned in C"
    } else {
      tumor_only <- sum(s$Tumor > 0 & s$Adjacent == 0)
      adjacent_only <- sum(s$Adjacent > 0 & s$Tumor == 0)
      note <- ""
    }
    single_cond_rows[[length(single_cond_rows) + 1]] <- data.frame(
      population = pop, cohort = co, tumor_only_donors = tumor_only,
      adjacent_only_donors = adjacent_only, note = note
    )
  }
}
single_cond <- do.call(rbind, single_cond_rows)
write.table(single_cond, file.path(out_dir, "single_condition_donors.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== Donors with cells in only one condition (candidates for unpaired/mixed model) ===\n")
print(single_cond, row.names = FALSE)

cat("\nDone. Tables written to results/ecm_population_choice/\n")
