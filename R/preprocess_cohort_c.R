## Individual pre-processing (integration step 1) for cohort C:
## NormalizeData + FindVariableFeatures on the symbol-indexed, post-QC-
## filtered checkpoint, computed independently of A/B.
##
## Input:  seurat_C_symbols.rds
## Output: seurat_C_normalized.rds
##         results/figures/integration/C/variable_features.png
##
## Run: Rscript R/preprocess_cohort_c.R

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

cat("Cohort C...\n")
obj_c <- readRDS(file.path(SCRATCH_DIR, "seurat_C_symbols.rds"))
obj_c <- preprocess_cohort(obj_c, "C", file.path(repo, "results", "figures", "integration", "C", "variable_features.png"))
saveRDS(obj_c, file.path(SCRATCH_DIR, "seurat_C_normalized.rds"))
cat("Saved seurat_C_normalized.rds\n")
