## Verification run (Section 3 of the submission task): reruns stages 03-06
## from the existing stage-01/02 checkpoints, in one clean Rscript session,
## recording wall-clock runtime and an approximate peak-memory figure per
## stage (gc()'s cumulative "max used" since the last gc(reset=TRUE), in Mb).
## Stages 01-02 are NOT rerun here (full reload/integration from raw inputs
## is hours of compute) - recorded as "not rerun after consolidation" in
## results/report_check.tsv instead, per the task's explicit fallback.
##
## This OVERWRITES each stage's own checkpoints/results files with freshly
## recomputed values (the existing scripts' normal behavior) - the resulting
## files are then compared against the report in a separate step
## (R/verify_report_numbers.R) to check true reproducibility, not merely
## that old files exist.
##
## Output: results/verification/runtime_memory_per_stage.tsv
## Run: Rscript R/verify_rerun_03_06.R

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

out_dir <- file.path(repo, "results", "verification")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stages <- list(
  "03-c-labeling" = c("annotate_c_by_cluster_majority.R", "verify_annotation_markers.R", "plot_abc_umap_by_population.R"),
  "04-ecm-overview-population-choice" = c("ecm_panel_by_population.R", "ecm_panel_by_population_per_cohort.R",
                                           "ecm_panel_by_population_per_condition.R", "ecm_population_choice_coverage.R"),
  "05-pseudobulk-models-robustness" = c("ecm_endothelial_pseudobulk_build.R", "ecm_endothelial_models.R",
                                         "ecm_endothelial_mesenchymal_flag.R", "ecm_endothelial_clean_models.R",
                                         "ecm_endothelial_add_Aonly_clean.R"),
  "06-linked-modalities" = c("multimodal_matched_subset.R", "multimodal_ecm_direction_check.R")
)

runtime_path <- file.path(out_dir, "runtime_memory_per_stage.tsv")
append_row <- function(df) {
  write.table(df, runtime_path, sep = "\t", row.names = FALSE, quote = FALSE,
              append = file.exists(runtime_path), col.names = !file.exists(runtime_path))
}
## Fresh file for this run (01/02 rows written immediately, so a crash in a
## later stage still leaves a complete record of everything that finished).
if (file.exists(runtime_path)) file.remove(runtime_path)
append_row(data.frame(stage = "01-load-qc", status = "not rerun after consolidation", runtime_sec = NA, peak_mem_mb = NA))
append_row(data.frame(stage = "02-integration", status = "not rerun after consolidation", runtime_sec = NA, peak_mem_mb = NA))

for (stage_name in names(stages)) {
  cat(sprintf("\n========== %s ==========\n", stage_name))
  gc(reset = TRUE, full = TRUE)
  t0 <- Sys.time()
  for (s in stages[[stage_name]]) {
    cat(sprintf("  source(R/%s)\n", s))
    source(file.path(repo, "R", s), local = new.env())
  }
  t1 <- Sys.time()
  gc_after <- gc()
  peak_mb <- sum(as.numeric(gc_after[, 6]))  # column 6 = "max used (Mb)", summed over Ncells + Vcells rows
  elapsed <- as.numeric(difftime(t1, t0, units = "secs"))
  cat(sprintf("  elapsed: %.1f sec, peak mem (approx, gc max used): %.1f Mb\n", elapsed, peak_mb))
  append_row(data.frame(stage = stage_name, status = "rerun (clean session, from existing checkpoints)",
                         runtime_sec = round(elapsed, 1), peak_mem_mb = round(peak_mb, 1)))
}

cat("\n=== Runtime / memory per stage ===\n")
print(read.delim(runtime_path), row.names = FALSE)
cat("\nDone.\n")
