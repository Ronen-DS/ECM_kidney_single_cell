## Driver: loads all 34 cohort A specimens (full matrices, restricted to
## author-annotated nuclei) and caches the result to disk, since this is the
## expensive step (~25 min - reading 34 raw droplet matrices).
## Run once: Rscript R/run_load_cohort_a.R

.args <- commandArgs(trailingOnly = FALSE)
.script_path <- sub("--file=", "", .args[grep("--file=", .args)])
repo <- normalizePath(file.path(dirname(.script_path), ".."))
readRenviron(file.path(repo, ".Renviron"))
source(file.path(repo, "config.R"))
source(file.path(repo, "R", "load_cohort_a.R"))

sample_manifest <- read.delim(SAMPLE_MANIFEST_PATH, stringsAsFactors = FALSE)

t0 <- Sys.time()
cohort_a <- load_cohort_a_all(sample_manifest)
cat("\nTotal elapsed:", as.numeric(Sys.time() - t0, units = "mins"), "min\n")

out_path <- file.path(SCRATCH_DIR, "cohort_a_loaded.rds")
saveRDS(cohort_a, out_path)
cat("Cached to:", out_path, "\n")

cat("\n=== Per-sample summary ===\n")
for (s in names(cohort_a)) {
  r <- cohort_a[[s]]
  cat(sprintf("%s  matched=%d  matrix=%dx%d\n", s, r$n_matched, nrow(r$matrix), ncol(r$matrix)))
}
cat("\nTotal cells across all samples:", sum(sapply(cohort_a, function(r) ncol(r$matrix))), "\n")
