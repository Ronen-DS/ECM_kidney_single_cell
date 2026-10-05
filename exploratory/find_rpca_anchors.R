## Classic Seurat anchor-based integration, step 2a: FindIntegrationAnchors
## (reduction="rpca") only, checkpointed separately from IntegrateData.
## Split out after IntegrateData crashed with a hard R-level allocation
## failure ('R_Calloc' could not allocate memory) on a first combined attempt -
## anchor-finding itself completed fine, so it shouldn't need to be redone if
## the (much heavier) IntegrateData step needs retrying with different
## parameters (see integrate_rpca_data.R).
##
## Input:  seurat_{A,B}_rpca_ready.rds
## Output: anchors_AB_rpca.rds
##
## Run: Rscript R/find_rpca_anchors.R

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

cat("Loading per-cohort RPCA-ready objects...\n")
obj_a <- readRDS(file.path(SCRATCH_DIR, "seurat_A_rpca_ready.rds"))
obj_b <- readRDS(file.path(SCRATCH_DIR, "seurat_B_rpca_ready.rds"))

cat("FindIntegrationAnchors(reduction = 'rpca', dims = 1:30)...\n")
anchors <- FindIntegrationAnchors(
  object.list = list(A = obj_a, B = obj_b),
  anchor.features = features_ab,
  reduction = "rpca",
  dims = 1:30,
  verbose = FALSE
)
rm(obj_a, obj_b)
gc(verbose = FALSE)

out_path <- file.path(SCRATCH_DIR, "anchors_AB_rpca.rds")
saveRDS(anchors, out_path)
cat("\nSaved:", out_path, "\n")
