## Re-derives seurat_C_symbols.rds from the updated seurat_C_filtered.rds
## (new non-MT-UMI-based QC thresholds - see filter_cohort_c.R) without
## touching A or B, which are unaffected by this change. Duplicate-symbol
## pairs themselves (24 for A, 34 for C) are a gene-level property of the
## reference, not the surviving cell set, so results/tables/
## duplicate_symbol_summary*.tsv from collapse_duplicate_symbols.R remain
## valid and are not regenerated here.
##
## Input:  seurat_C_filtered.rds
## Output: seurat_C_symbols.rds
##
## Run: Rscript R/recollapse_cohort_c.R

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
source(file.path(repo, "R", "symbol_utils.R"))
library(Seurat)
library(Matrix)

cat("Cohort C...\n")
obj_c <- readRDS(file.path(SCRATCH_DIR, "seurat_C_filtered.rds"))
features_c <- obj_c@misc$features
counts_c <- GetAssayData(obj_c, assay = "RNA", layer = "counts")
new_counts_c <- collapse_symbol_duplicates(counts_c, features_c)

meta <- obj_c@meta.data
meta$nCount_RNA <- NULL
meta$nFeature_RNA <- NULL
obj_c_symbols <- CreateSeuratObject(counts = new_counts_c, meta.data = meta, project = "cohort_C")
obj_c_symbols@misc$cohort <- "C"
obj_c_symbols@misc$features_pre_collapse <- features_c

saveRDS(obj_c_symbols, file.path(SCRATCH_DIR, "seurat_C_symbols.rds"))
cat(sprintf("Saved seurat_C_symbols.rds: %d genes x %d cells\n", nrow(obj_c_symbols), ncol(obj_c_symbols)))
