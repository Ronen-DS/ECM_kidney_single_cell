## Prepares the unintegrated A+B+C baseline. Regenerates
## SelectIntegrationFeatures fresh for A+B+C (cohort C's population changed
## after re-filtering, so the stale features file from the original,
## abandoned 3-cohort attempt is not reused), then subsets EACH cohort to
## just those 2000 features BEFORE merging - not after, unlike the original
## attempt that caused a near-OOM crash by merging full ~45k-gene objects
## first. This should make peak memory far lower even with 3 cohorts loaded
## simultaneously.
##
## Input:  seurat_{A,B,C}_normalized.rds
## Output: merged_ABC.rds
##         results/tables/integration/integration_features_ABC.txt
##
## Run: Rscript R/prepare_integration_abc.R

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
source(file.path(repo, "R", "integration_utils.R"))
library(Seurat)

cat("Loading normalized per-cohort checkpoints (A, B, C)...\n")
obj_a <- readRDS(file.path(SCRATCH_DIR, "seurat_A_normalized.rds"))
obj_b <- readRDS(file.path(SCRATCH_DIR, "seurat_B_normalized.rds"))
obj_c <- readRDS(file.path(SCRATCH_DIR, "seurat_C_normalized.rds"))
cat(sprintf("  A: %d x %d | B: %d x %d | C: %d x %d\n",
            nrow(obj_a), ncol(obj_a), nrow(obj_b), ncol(obj_b), nrow(obj_c), ncol(obj_c)))

cat("SelectIntegrationFeatures([A,B,C])...\n")
features_abc <- SelectIntegrationFeatures(object.list = list(A = obj_a, B = obj_b, C = obj_c), nfeatures = 2000)

tables_dir <- file.path(repo, "results", "tables", "integration")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
writeLines(features_abc, file.path(tables_dir, "integration_features_ABC.txt"))
cat(sprintf("  %d integration features\n", length(features_abc)))

cat("Subsetting each cohort to the shared feature set BEFORE merging...\n")
obj_a <- subset(obj_a, features = features_abc)
obj_b <- subset(obj_b, features = features_abc)
obj_c <- subset(obj_c, features = features_abc)
cat(sprintf("  Subsetted: A %d x %d | B %d x %d | C %d x %d\n",
            nrow(obj_a), ncol(obj_a), nrow(obj_b), ncol(obj_b), nrow(obj_c), ncol(obj_c)))

cat("Merging [A+B+C]...\n")
merged_abc <- merge_for_integration(list(A = obj_a, B = obj_b, C = obj_c))
rm(obj_a, obj_b, obj_c)
gc(verbose = FALSE)
cat(sprintf("  Merged: %d genes x %d cells\n", nrow(merged_abc), ncol(merged_abc)))
cat("  Layers:", paste(Layers(merged_abc[["RNA"]]), collapse = ", "), "\n")

out_path <- file.path(SCRATCH_DIR, "merged_ABC.rds")
saveRDS(merged_abc, out_path)
cat("\nSaved:", out_path, "\n")
