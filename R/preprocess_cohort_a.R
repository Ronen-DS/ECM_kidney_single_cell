## Individual pre-processing (integration step 1) for cohort A:
## NormalizeData + FindVariableFeatures on the symbol-indexed, post-QC-
## filtered checkpoint, computed independently of B/C.
##
## Input:  seurat_A_symbols.rds
## Output: seurat_A_normalized.rds
##         results/figures/integration/A/variable_features.png
##
## Run: Rscript R/preprocess_cohort_a.R

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

cat("Cohort A...\n")
obj_a <- readRDS(file.path(SCRATCH_DIR, "seurat_A_symbols.rds"))
obj_a <- preprocess_cohort(obj_a, "A", file.path(repo, "results", "figures", "integration", "A", "variable_features.png"))
saveRDS(obj_a, file.path(SCRATCH_DIR, "seurat_A_normalized.rds"))
cat("Saved seurat_A_normalized.rds\n")
