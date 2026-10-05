## Individual pre-processing (integration step 1) for cohort B:
## NormalizeData + FindVariableFeatures on the symbol-indexed, post-QC-
## filtered checkpoint, computed independently of A/C.
##
## Input:  seurat_B_symbols.rds
## Output: seurat_B_normalized.rds
##         results/figures/integration/B/variable_features.png
##
## Run: Rscript R/preprocess_cohort_b.R

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

cat("Cohort B...\n")
obj_b <- readRDS(file.path(SCRATCH_DIR, "seurat_B_symbols.rds"))
obj_b <- preprocess_cohort(obj_b, "B", file.path(repo, "results", "figures", "integration", "B", "variable_features.png"))
saveRDS(obj_b, file.path(SCRATCH_DIR, "seurat_B_normalized.rds"))
cat("Saved seurat_B_normalized.rds\n")
