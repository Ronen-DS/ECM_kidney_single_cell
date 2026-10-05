## Combine and prepare for integration (integration step 2): SelectIntegrationFeatures()
## for [A+B], then merge cohorts A and B into one Seurat object. Seurat v5's
## merge() keeps each input object's counts/data as separate per-cohort
## layers automatically (see integration_utils.R), so a later
## IntegrateLayers() call already has batches to work with. No batch
## correction happens in this script - the merged object IS the unintegrated
## baseline that an integration method will later be compared against.
##
## Cohort C is deliberately NOT loaded here (dropped from scope - A+B is
## sufficient) so this script's peak memory footprint is just A+B, not
## A+B+C - loading all three normalized checkpoints upfront (each already
## holding both counts and log-normalized data layers for 120-160k cells)
## was the actual cause of a near-OOM run of an earlier 3-cohort version of
## this script.
##
## Input:  seurat_{A,B}_normalized.rds
## Output: merged_AB.rds
##         results/tables/integration/integration_features_AB.txt
##         results/tables/integration/integration_prep_summary.tsv
##
## Run: Rscript R/prepare_integration.R

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

ecm_genes <- read.delim(ECM_PANEL_PATH, stringsAsFactors = FALSE)$gene_symbol

cat("Loading normalized per-cohort checkpoints (A, B only)...\n")
obj_a <- readRDS(file.path(SCRATCH_DIR, "seurat_A_normalized.rds"))
obj_b <- readRDS(file.path(SCRATCH_DIR, "seurat_B_normalized.rds"))

tables_dir <- file.path(repo, "results", "tables", "integration")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

## ---- [A+B] -------------------------------------------------------------
cat("\n[A+B] SelectIntegrationFeatures...\n")
ab_list <- list(A = obj_a, B = obj_b)
ab_summary <- integration_feature_summary(ab_list, ecm_genes)
writeLines(ab_summary$features, file.path(tables_dir, "integration_features_AB.txt"))
cat(sprintf("  %d integration features; %d/%d ECM panel genes included (%s)\n",
            ab_summary$n_features, ab_summary$n_ecm_included, length(ecm_genes),
            paste(ab_summary$ecm_genes_included, collapse = ", ")))

n_genes_union <- length(unique(c(rownames(obj_a), rownames(obj_b))))
n_cells <- ncol(obj_a) + ncol(obj_b)
rm(ab_list)

cat("Merging [A+B]...\n")
merged_ab <- merge_for_integration(list(A = obj_a, B = obj_b))
rm(obj_a, obj_b)
gc(verbose = FALSE)
cat("  Layers:", paste(Layers(merged_ab[["RNA"]]), collapse = ", "), "\n")
saveRDS(merged_ab, file.path(SCRATCH_DIR, "merged_AB.rds"))
cat(sprintf("  Saved merged_AB.rds: %d genes x %d cells\n", nrow(merged_ab), ncol(merged_ab)))

## ---- Summary table ---------------------------------------------------------
summary_tbl <- data.frame(
  combination = "A+B",
  n_integration_features = ab_summary$n_features,
  n_ecm_genes_included = ab_summary$n_ecm_included,
  n_cells = n_cells,
  n_genes_union = n_genes_union
)
write.table(summary_tbl, file.path(tables_dir, "integration_prep_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n")
print(summary_tbl, row.names = FALSE)
