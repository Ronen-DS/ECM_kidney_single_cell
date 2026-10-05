## Verification (Section 3, correction): reruns ONLY stage 05 (pseudobulk
## models + robustness) as its own fresh, short-lived R process - not part of
## a long multi-stage session - to get an accurate runtime/peak-memory
## reading. The previous stage-05 timing (23656.2 sec / 6.6h) was an artifact
## of a long gap during a different run's manual polling, not real compute
## time, and needed a clean, uninterrupted measurement.
##
## Uses the existing merged_ABC_annotated.rds checkpoint (stage 03's output,
## itself unaffected - the stage 01/02 raw-input rerun attempted after this
## stage-05 reading was last taken crashed with std::bad_alloc during
## IntegrateData and was abandoned per the user's decision; stage 02's
## checkpoint was never overwritten, so merged_ABC_annotated.rds on disk is
## still the same, already-verified one).
##
## Updates results/verification/runtime_memory_per_stage.tsv: corrects the
## 05-pseudobulk-models-robustness row, and adds an honest 02-integration row
## documenting the crash (rather than leaving that stage's row silently
## absent).
##
## Run: Rscript R/verify_rerun_05_only.R

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
runtime_path <- file.path(out_dir, "runtime_memory_per_stage.tsv")

stage_scripts <- c(
  "ecm_endothelial_pseudobulk_build.R",
  "ecm_endothelial_models.R",
  "ecm_endothelial_mesenchymal_flag.R",
  "ecm_endothelial_clean_models.R",
  "ecm_endothelial_add_Aonly_clean.R"
)

cat("========== 05-pseudobulk-models-robustness (fresh process, clean timing) ==========\n")
gc(reset = TRUE, full = TRUE)
t0 <- Sys.time()
for (s in stage_scripts) {
  cat(sprintf("  [%s] source(R/%s)\n", format(Sys.time(), "%H:%M:%S"), s))
  source(file.path(repo, "R", s), local = new.env())
}
t1 <- Sys.time()
gc_after <- gc()
peak_mb <- sum(as.numeric(gc_after[, 6]))
elapsed <- as.numeric(difftime(t1, t0, units = "secs"))
cat(sprintf("\n[05-pseudobulk-models-robustness] elapsed: %.1f sec (%.1f min), peak mem (approx, gc max used): %.1f Mb\n",
            elapsed, elapsed / 60, peak_mb))

## --- Update runtime_memory_per_stage.tsv: fix 05's row, add an honest 02 row ---
existing <- read.delim(runtime_path, stringsAsFactors = FALSE)
existing <- existing[existing$stage != "05-pseudobulk-models-robustness", ]
new_05 <- data.frame(stage = "05-pseudobulk-models-robustness",
                      status = "rerun (fresh standalone process, clean uninterrupted timing)",
                      runtime_sec = round(elapsed, 1), peak_mem_mb = round(peak_mb, 1))
combined <- rbind(existing, new_05)

if (!("02-integration" %in% combined$stage)) {
  combined <- rbind(combined, data.frame(
    stage = "02-integration",
    status = "attempted from raw inputs, CRASHED (std::bad_alloc during IntegrateData, k.weight=30) after ~5h in one long session - not completed; abandoned, stage 02 checkpoint unchanged from the earlier successful run",
    runtime_sec = NA, peak_mem_mb = NA))
}

order_ref <- c("01-load-qc", "02-integration", "03-c-labeling", "04-ecm-overview-population-choice",
               "05-pseudobulk-models-robustness", "06-linked-modalities")
combined <- combined[match(order_ref, combined$stage), ]
combined <- combined[!is.na(combined$stage), ]
write.table(combined, runtime_path, sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== Updated runtime_memory_per_stage.tsv ===\n")
print(combined, row.names = FALSE)
cat("\nDone.\n")
