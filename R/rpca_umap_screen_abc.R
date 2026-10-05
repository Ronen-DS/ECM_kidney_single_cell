## Classic Seurat anchor-based integration for A+B+C, final step: on the
## "integrated" assay, ScaleData + RunPCA + RunUMAP, straight to a by-cohort
## DimPlot - this is figure 2 of the user's 2-figure comparison (figure 1 is
## the unintegrated A+B+C baseline from baseline_abc_umap.R).
##
## Input:  integrated_ABC_rpca.rds
## Output: merged_ABC_rpca_umap.rds
##         results/figures/integration/ABC_umap_rpca_by_cohort.png
##
## Run: Rscript R/rpca_umap_screen_abc.R

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
library(ggplot2)

n_dims <- 10

cat("Loading integrated_ABC_rpca.rds...\n")
integrated <- readRDS(file.path(SCRATCH_DIR, "integrated_ABC_rpca.rds"))
DefaultAssay(integrated) <- "integrated"

cat("ScaleData + RunPCA on integrated assay...\n")
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)

cat(sprintf("RunUMAP(dims = 1:%d)...\n", n_dims))
integrated <- RunUMAP(integrated, dims = 1:n_dims, verbose = FALSE)

fig_dir <- file.path(repo, "results", "figures", "integration")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
abc_colors <- c(A = "#F8766D", B = "#00BFC4", C = "#7CAE00")
p <- DimPlot(integrated, reduction = "umap", group.by = "cohort") +
  scale_color_manual(values = abc_colors) +
  ggtitle("Integrated [A+B+C] (RPCA anchors): UMAP by cohort")
ggsave(file.path(fig_dir, "ABC_umap_rpca_by_cohort.png"), p, width = 8, height = 6, dpi = 150, bg = "white")

out_path <- file.path(SCRATCH_DIR, "merged_ABC_rpca_umap.rds")
saveRDS(integrated, out_path)
cat("\nSaved:", out_path, "\n")
cat("Wrote results/figures/integration/ABC_umap_rpca_by_cohort.png\n")
