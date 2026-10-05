## Qualitative integration check: re-colors the already-computed UMAPs
## (unintegrated, Harmony, RPCA) by a harmonized cell-population label
## (metadata/label_mapping.tsv), instead of by cohort. If integration
## worked, cells sharing a population label should co-localize regardless of
## cohort; if a method just merges everything indiscriminately, population
## islands would look scrambled instead of clean. Qualitative only - no
## metrics computed (see integration_utils.R::add_major_population() for the
## join logic, including B's cross_population_conflict handling).
##
## Input:  merged_AB_umap.rds, merged_AB_harmony_umap.rds, merged_AB_rpca_final.rds
##         metadata/label_mapping.tsv
## Output: results/figures/integration/AB_umap_unintegrated_by_population.png
##         results/figures/integration/AB_umap_harmony_by_population.png
##         results/figures/integration/AB_umap_rpca_by_population.png
##
## Run: Rscript R/plot_by_population.R

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

label_mapping <- read.delim(LABEL_MAPPING_PATH, stringsAsFactors = FALSE)
fig_dir <- file.path(repo, "results", "figures", "integration")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## A fixed color per population, shared across all 3 plots, so the same
## population is the same color in every panel for a fair visual comparison.
populations <- sort(unique(label_mapping$population))
pop_colors <- setNames(scales::hue_pal()(length(populations)), populations)
pop_colors["Unresolved"] <- "grey80"

plot_one <- function(rds_name, reduction, title, out_name) {
  cat(sprintf("\n%s (reduction = %s)...\n", rds_name, reduction))
  obj <- readRDS(file.path(SCRATCH_DIR, rds_name))
  obj <- add_major_population(obj, label_mapping)
  p <- DimPlot(obj, reduction = reduction, group.by = "major_population") +
    scale_color_manual(values = pop_colors) +
    ggtitle(title)
  ggsave(file.path(fig_dir, out_name), p, width = 8, height = 6, dpi = 150, bg = "white")
  cat("  Wrote", out_name, "\n")
  rm(obj)
  gc(verbose = FALSE)
}

plot_one("merged_AB_umap.rds", "umap",
         "Unintegrated [A+B]: UMAP by harmonized cell population",
         "AB_umap_unintegrated_by_population.png")

plot_one("merged_AB_harmony_umap.rds", "umap.harmony",
         "Integrated [A+B] (Harmony): UMAP by harmonized cell population",
         "AB_umap_harmony_by_population.png")

plot_one("merged_AB_rpca_final.rds", "umap",
         "Integrated [A+B] (RPCA anchors): UMAP by harmonized cell population",
         "AB_umap_rpca_by_population.png")

cat("\nDone.\n")
