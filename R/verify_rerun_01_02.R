## Verification run (Section 3, continued): reruns stages 01-02 from raw
## inputs (not from any existing checkpoint - none of these scripts have
## internal skip-if-exists logic, confirmed before launching, so sourcing
## them always reloads/recomputes from the true source files under
## DATA_ROOT). Updates results/verification/runtime_memory_per_stage.tsv,
## replacing the "not rerun after consolidation" placeholder rows for
## 01-load-qc and 02-integration with real timing/memory, while leaving the
## already-recorded 03-06 rows untouched.
##
## Run: Rscript R/verify_rerun_01_02.R

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
runtime_path <- file.path(out_dir, "runtime_memory_per_stage.tsv")

stages <- list(
  "01-load-qc" = c(
    "run_load_cohort_a.R", "build_seurat_a.R",
    "build_seurat_b.R",
    "build_seurat_c.R",
    "qc_cohort_a.R", "qc_cohort_b.R", "qc_cohort_c.R",
    "derive_c_mt_threshold.R",
    "filter_cohort_a.R", "filter_cohort_b.R", "filter_cohort_c.R",
    "collapse_duplicate_symbols.R",
    "recollapse_cohort_c.R",
    "preprocess_cohort_a.R", "preprocess_cohort_b.R", "preprocess_cohort_c.R",
    "qc_population_retention.R", "qc_population_retention_A5pct.R",
    "qc_combined_ABC_figure.R"
  ),
  "02-integration" = c(
    "prepare_integration_abc.R",
    "baseline_abc_umap.R",
    "prepare_rpca_pca_abc.R",
    "find_rpca_anchors_abc.R",
    "integrate_rpca_data_abc.R",
    "rpca_pca_only_abc.R",
    "rpca_umap_only_abc.R",
    "rpca_clusters_umap_abc.R",
    "rpca_cluster_res03_abc.R",
    "crosstab_res03_population_abc.R"
  )
)

## Preserve existing 03-06 rows, drop the old 01/02 placeholder rows.
existing <- if (file.exists(runtime_path)) read.delim(runtime_path, stringsAsFactors = FALSE) else NULL
keep_existing <- if (!is.null(existing)) existing[!(existing$stage %in% names(stages)), ] else NULL

new_rows <- list()
for (stage_name in names(stages)) {
  cat(sprintf("\n========== %s (from raw inputs) ==========\n", stage_name))
  gc(reset = TRUE, full = TRUE)
  t0 <- Sys.time()
  for (s in stages[[stage_name]]) {
    cat(sprintf("  [%s] source(R/%s)\n", format(Sys.time(), "%H:%M:%S"), s))
    source(file.path(repo, "R", s), local = new.env())
  }
  t1 <- Sys.time()
  gc_after <- gc()
  peak_mb <- sum(as.numeric(gc_after[, 6]))
  elapsed <- as.numeric(difftime(t1, t0, units = "secs"))
  cat(sprintf("  [%s] elapsed: %.1f sec (%.1f min), peak mem (approx, gc max used): %.1f Mb\n",
              stage_name, elapsed, elapsed / 60, peak_mb))
  new_rows[[stage_name]] <- data.frame(stage = stage_name, status = "rerun from raw inputs (clean session)",
                                        runtime_sec = round(elapsed, 1), peak_mem_mb = round(peak_mb, 1))

  ## Write incrementally after EACH stage, so a crash during 02 still leaves
  ## 01's hard-won timing on disk.
  combined <- rbind(keep_existing, do.call(rbind, new_rows))
  order_ref <- c("01-load-qc", "02-integration", "03-c-labeling", "04-ecm-overview-population-choice",
                 "05-pseudobulk-models-robustness", "06-linked-modalities")
  combined <- combined[match(order_ref, combined$stage), ]
  combined <- combined[!is.na(combined$stage), ]
  write.table(combined, runtime_path, sep = "\t", row.names = FALSE, quote = FALSE)
}

cat("\n=== Runtime / memory per stage (final) ===\n")
print(read.delim(runtime_path), row.names = FALSE)
cat("\nDone.\n")
