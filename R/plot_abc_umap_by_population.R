## Re-colors the existing joint A+B+C RPCA-integrated UMAP (already computed
## in rpca_umap_only_abc.R) by harmonized cell population (label_mapping.tsv)
## for cohorts A and B - the ones with real annotations - while cohort C's
## cells keep their own distinct cohort color instead of falling into a
## generic "Unresolved" bucket, since C has no supplied cell-type
## annotations at all (not a case of ambiguous/conflicting labels like some
## A/B cells - there's simply nothing to map).
##
## Input:  merged_ABC_rpca_umap.rds, metadata/label_mapping.tsv
## Output: results/figures/integration/ABC_umap_rpca_by_population.png
##
## Run: Rscript R/plot_abc_umap_by_population.R

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
library(Seurat)
library(ggplot2)

label_mapping <- read.delim(LABEL_MAPPING_PATH, stringsAsFactors = FALSE)

cat("Loading merged_ABC_rpca_umap.rds...\n")
obj <- readRDS(file.path(SCRATCH_DIR, "merged_ABC_rpca_umap.rds"))
obj <- add_major_population(obj, label_mapping) # C rows fall to "Unresolved" here (no A/B labels apply) - overridden below

is_c <- obj$cohort == "C"
plot_group <- obj$major_population
plot_group[is_c] <- "Cohort C"
obj$plot_group <- plot_group
cat(sprintf("Cohort C cells kept in their own category: %d\n", sum(is_c)))

populations <- sort(unique(label_mapping$population))
pop_colors <- setNames(scales::hue_pal()(length(populations)), populations)
pop_colors["Unresolved"] <- "grey80"
## NOT the same olive-green used for cohort C elsewhere - that collided
## visually with Malignant_epithelial's hue-palette color here, defeating
## the point of telling them apart. Black is distinct from every population
## color and from Unresolved's light grey.
pop_colors["Cohort C"] <- "black"

fig_dir <- file.path(repo, "results", "figures", "integration")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
p <- DimPlot(obj, reduction = "umap", group.by = "plot_group") +
  scale_color_manual(values = pop_colors) +
  ggtitle("Integrated [A+B+C] (RPCA anchors): UMAP by cell population (A/B) / cohort C")
ggsave(file.path(fig_dir, "ABC_umap_rpca_by_population.png"), p, width = 8.5, height = 6, dpi = 150, bg = "white")

cat("Wrote results/figures/integration/ABC_umap_rpca_by_population.png\n")
