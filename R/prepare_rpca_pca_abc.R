## Classic Seurat anchor-based integration for A+B+C, step 1: per-cohort
## ScaleData + RunPCA, computed INDIVIDUALLY on A, B, and C, using the fresh
## A+B+C SelectIntegrationFeatures() (integration_features_ABC.txt - C's
## population changed since re-filtering, so this is NOT the A+B-only
## feature list reused).
##
## Input:  seurat_{A,B,C}_normalized.rds, results/tables/integration/integration_features_ABC.txt
## Output: seurat_A_rpca_ready_abc.rds, seurat_B_rpca_ready_abc.rds, seurat_C_rpca_ready_abc.rds
##
## Run: Rscript R/prepare_rpca_pca_abc.R

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

features_abc <- readLines(file.path(repo, "results", "tables", "integration", "integration_features_ABC.txt"))
cat(sprintf("Using %d integration features.\n", length(features_abc)))

for (cohort in c("A", "B", "C")) {
  cat(sprintf("\nCohort %s...\n", cohort))
  in_path <- file.path(SCRATCH_DIR, sprintf("seurat_%s_normalized.rds", cohort))
  obj <- readRDS(in_path)
  obj <- subset(obj, features = features_abc)
  cat(sprintf("  Subsetted: %d genes x %d cells\n", nrow(obj), ncol(obj)))
  obj <- ScaleData(obj, features = features_abc, verbose = FALSE)
  obj <- RunPCA(obj, features = features_abc, npcs = 30, verbose = FALSE)
  out_path <- file.path(SCRATCH_DIR, sprintf("seurat_%s_rpca_ready_abc.rds", cohort))
  saveRDS(obj, out_path)
  cat("  Saved:", out_path, "\n")
  rm(obj)
  gc(verbose = FALSE)
}

cat("\nDone.\n")
