## Classic Seurat anchor-based integration for A+B+C, step 2b: IntegrateData.
## Uses k.weight = 30 (not the default 100) from the start - the A+B version
## of this exact step crashed with a genuine R-level 'R_Calloc could not
## allocate memory' error at the default, and k.weight=30 was the fix that
## worked there. No reason to repeat that guaranteed-to-fail first attempt.
##
## Input:  anchors_ABC_rpca.rds
## Output: integrated_ABC_rpca.rds
##
## Run: Rscript R/integrate_rpca_data_abc.R

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

K_WEIGHT <- 30

cat("Loading anchors_ABC_rpca.rds...\n")
anchors <- readRDS(file.path(SCRATCH_DIR, "anchors_ABC_rpca.rds"))

cat(sprintf("IntegrateData(dims = 1:30, k.weight = %d)...\n", K_WEIGHT))
integrated <- IntegrateData(anchorset = anchors, dims = 1:30, k.weight = K_WEIGHT, verbose = FALSE)
rm(anchors)
gc(verbose = FALSE)

out_path <- file.path(SCRATCH_DIR, "integrated_ABC_rpca.rds")
saveRDS(integrated, out_path)
cat("\nSaved:", out_path, "\n")
cat(sprintf("Integrated assay: %d features x %d cells\n", nrow(integrated), ncol(integrated)))
