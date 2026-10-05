## Re-saves the [A+B] elbow plot with an explicit white background - the
## first save used ggsave()'s default transparent background, which some
## viewers render as solid black instead of transparent/white. Loads the
## already-computed merged_AB_pca.rds checkpoint (cheap) rather than
## re-running ScaleData/RunPCA.
##
## Run: Rscript R/replot_elbow_ab.R

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

merged_ab <- readRDS(file.path(SCRATCH_DIR, "merged_AB_pca.rds"))

fig_dir <- file.path(repo, "results", "figures", "integration")
p <- ElbowPlot(merged_ab, ndims = 50) + ggtitle("Unintegrated [A+B]: PCA elbow plot")
ggsave(file.path(fig_dir, "AB_elbow_plot.png"), p, width = 7, height = 5, dpi = 150,
       bg = "white")
cat("Re-saved results/figures/integration/AB_elbow_plot.png with a white background.\n")
