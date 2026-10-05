## Classic Seurat anchor-based integration, step 1: per-cohort ScaleData +
## RunPCA, computed INDIVIDUALLY on A and B (not jointly) - FindIntegrationAnchors
## with reduction="rpca" requires each dataset's own PCA embedding to project
## the other dataset into for anchor search. Uses the same 2000
## SelectIntegrationFeatures() genes already chosen in prepare_integration.R,
## for consistency with the (unsuccessful) Harmony attempt.
##
## Input:  seurat_{A,B}_normalized.rds, results/tables/integration/integration_features_AB.txt
## Output: seurat_A_rpca_ready.rds, seurat_B_rpca_ready.rds
##
## Run: Rscript R/prepare_rpca_pca.R

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

features_ab <- readLines(file.path(repo, "results", "tables", "integration", "integration_features_AB.txt"))
cat(sprintf("Using %d integration features.\n", length(features_ab)))

for (cohort in c("A", "B")) {
  cat(sprintf("\nCohort %s...\n", cohort))
  in_path <- file.path(SCRATCH_DIR, sprintf("seurat_%s_normalized.rds", cohort))
  obj <- readRDS(in_path)
  obj <- subset(obj, features = features_ab)
  cat(sprintf("  Subsetted: %d genes x %d cells\n", nrow(obj), ncol(obj)))
  obj <- ScaleData(obj, features = features_ab, verbose = FALSE)
  obj <- RunPCA(obj, features = features_ab, npcs = 30, verbose = FALSE)
  out_path <- file.path(SCRATCH_DIR, sprintf("seurat_%s_rpca_ready.rds", cohort))
  saveRDS(obj, out_path)
  cat("  Saved:", out_path, "\n")
  rm(obj)
  gc(verbose = FALSE)
}

cat("\nDone.\n")
