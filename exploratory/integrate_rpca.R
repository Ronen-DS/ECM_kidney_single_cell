## Classic Seurat anchor-based integration, step 2: FindIntegrationAnchors
## (reduction="rpca" - faster/lighter than full CCA, chosen given this
## machine's repeated near-OOM issues at 245k total cells) + IntegrateData,
## producing a new "integrated" assay holding batch-corrected values for the
## 2000 anchor features across all A+B cells.
##
## Input:  seurat_{A,B}_rpca_ready.rds
## Output: integrated_AB_rpca.rds
##
## Run: Rscript R/integrate_rpca.R

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

cat("IntegrateData(dims = 1:30)...\n")
integrated <- IntegrateData(anchorset = anchors, dims = 1:30, verbose = FALSE)
rm(anchors)
gc(verbose = FALSE)

out_path <- file.path(SCRATCH_DIR, "integrated_AB_rpca.rds")
saveRDS(integrated, out_path)
cat("\nSaved:", out_path, "\n")
cat(sprintf("Integrated assay: %d features x %d cells\n", nrow(integrated), ncol(integrated)))
