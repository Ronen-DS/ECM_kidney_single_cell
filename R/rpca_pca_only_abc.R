## Split out from rpca_umap_screen_abc.R after two attempts at the combined
## script were killed by the harness's background time limit (not a crash or
## memory issue) following a session restart - isolating ScaleData+RunPCA
## (historically fast, a few minutes) from RunUMAP (the long step) to
## checkpoint progress and narrow down whether the time-limit issue is
## specific to sustained long-running background tasks post-restart.
##
## Input:  integrated_ABC_rpca.rds
## Output: integrated_ABC_rpca_pca.rds
##         results/figures/integration/ABC_elbow_plot_rpca.png
##
## Run: Rscript R/rpca_pca_only_abc.R

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

cat("Loading integrated_ABC_rpca.rds...\n")
integrated <- readRDS(file.path(SCRATCH_DIR, "integrated_ABC_rpca.rds"))
DefaultAssay(integrated) <- "integrated"

cat("ScaleData + RunPCA on integrated assay...\n")
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)

fig_dir <- file.path(repo, "results", "figures", "integration")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
p_elbow <- ElbowPlot(integrated, ndims = 30) + ggtitle("Integrated [A+B+C] (RPCA anchors): elbow plot")
ggsave(file.path(fig_dir, "ABC_elbow_plot_rpca.png"), p_elbow, width = 7, height = 5, dpi = 150, bg = "white")

out_path <- file.path(SCRATCH_DIR, "integrated_ABC_rpca_pca.rds")
saveRDS(integrated, out_path)
cat("\nSaved:", out_path, "\n")
