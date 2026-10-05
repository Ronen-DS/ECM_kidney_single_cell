## Classic Seurat anchor-based integration for A+B+C, step 2a:
## FindIntegrationAnchors(reduction="rpca") across all three cohorts
## (computes anchors for all 3 pairs: A-B, A-C, B-C), checkpointed
## separately from IntegrateData - mirrors the A+B pipeline's split, done
## from the start this time since IntegrateData is the step known to be
## memory-risky, not anchor-finding.
##
## Input:  seurat_{A,B,C}_rpca_ready_abc.rds
## Output: anchors_ABC_rpca.rds
##
## Run: Rscript R/find_rpca_anchors_abc.R

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

cat("Loading per-cohort RPCA-ready objects (A, B, C)...\n")
obj_a <- readRDS(file.path(SCRATCH_DIR, "seurat_A_rpca_ready_abc.rds"))
obj_b <- readRDS(file.path(SCRATCH_DIR, "seurat_B_rpca_ready_abc.rds"))
obj_c <- readRDS(file.path(SCRATCH_DIR, "seurat_C_rpca_ready_abc.rds"))

cat("FindIntegrationAnchors(reduction = 'rpca', dims = 1:30, across A/B/C)...\n")
anchors <- FindIntegrationAnchors(
  object.list = list(A = obj_a, B = obj_b, C = obj_c),
  anchor.features = features_abc,
  reduction = "rpca",
  dims = 1:30,
  verbose = FALSE
)
rm(obj_a, obj_b, obj_c)
gc(verbose = FALSE)

out_path <- file.path(SCRATCH_DIR, "anchors_ABC_rpca.rds")
saveRDS(anchors, out_path)
cat("\nSaved:", out_path, "\n")
