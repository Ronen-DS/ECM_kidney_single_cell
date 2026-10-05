## Classic Seurat anchor-based integration, step 2b: IntegrateData, split out
## from anchor-finding after a first combined attempt crashed with a hard
## R-level allocation failure ('R_Calloc' could not allocate memory) inside
## IntegrateData's weighted-correction step. Uses a reduced k.weight (30,
## down from the default 100) - Seurat's own guidance for large datasets is
## to reduce k.weight to shrink the per-cell weighting-neighborhood matrices
## IntegrateData builds internally, which is what actually exhausted memory
## last time (anchor-finding itself completed without issue).
##
## Input:  anchors_AB_rpca.rds
## Output: integrated_AB_rpca.rds
##
## Run: Rscript R/integrate_rpca_data.R

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

K_WEIGHT <- 30 # default is 100; reduced to shrink IntegrateData's internal weighting matrices

cat("Loading anchors_AB_rpca.rds...\n")
anchors <- readRDS(file.path(SCRATCH_DIR, "anchors_AB_rpca.rds"))

cat(sprintf("IntegrateData(dims = 1:30, k.weight = %d)...\n", K_WEIGHT))
integrated <- IntegrateData(anchorset = anchors, dims = 1:30, k.weight = K_WEIGHT, verbose = FALSE)
rm(anchors)
gc(verbose = FALSE)

out_path <- file.path(SCRATCH_DIR, "integrated_AB_rpca.rds")
saveRDS(integrated, out_path)
cat("\nSaved:", out_path, "\n")
cat(sprintf("Integrated assay: %d features x %d cells\n", nrow(integrated), ncol(integrated)))
