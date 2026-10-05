## Helpers for per-cohort pre-processing and the [A+B]/[A+B+C] "prepare for
## integration" step. Kept notebook-compatible: pure functions taking
## explicit arguments (object, label, output path) rather than reaching for
## commandArgs-derived paths, so they drop into a notebook cell unchanged.

library(Seurat)
library(ggplot2)

#' NormalizeData + FindVariableFeatures on ONE cohort's Seurat object,
#' independently (no cross-cohort information used here) - the "individual
#' pre-processing" step that SelectIntegrationFeatures() later needs each
#' object in its list to already have. Writes a VariableFeaturePlot with the
#' top 10 HVGs labeled. Returns the updated object; the caller decides
#' whether/where to checkpoint it.
preprocess_cohort <- function(obj, cohort_label, fig_path, nfeatures = 2000) {
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)

  top10 <- head(VariableFeatures(obj), 10)
  plot1 <- VariableFeaturePlot(obj)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  combined <- (plot1 + plot2) +
    patchwork::plot_annotation(title = sprintf("Cohort %s: highly variable features", cohort_label))
  dir.create(dirname(fig_path), showWarnings = FALSE, recursive = TRUE)
  ggsave(fig_path, combined, width = 12, height = 5, dpi = 150)

  cat(sprintf("  Cohort %s top 10 HVGs: %s\n", cohort_label, paste(top10, collapse = ", ")))
  obj
}

#' SelectIntegrationFeatures() for a list of already-preprocessed cohort
#' objects, plus how many of a reference gene panel (e.g. the ECM panel) made
#' the cut - a sanity signal, not a requirement (integration features are
#' chosen for overall variance structure, not to preserve any one panel).
integration_feature_summary <- function(object_list, ecm_genes, nfeatures = 2000) {
  features <- SelectIntegrationFeatures(object.list = object_list, nfeatures = nfeatures)
  list(
    features = features,
    n_features = length(features),
    ecm_genes_included = intersect(ecm_genes, features),
    n_ecm_included = length(intersect(ecm_genes, features))
  )
}

#' Merges a named list of cohort Seurat objects into one. Seurat v5's
#' merge() already keeps each input object's counts/data as SEPARATE layers
#' by default (e.g. "counts.cohort_A", "counts.cohort_B", one pair per input
#' object) rather than flattening them into a single shared layer - this IS
#' the "layer-splitting by specimen/batch" the assignment asks for, done
#' automatically by merge() itself, so no extra split() call is needed (and
#' calling one errors, since the layers are already split). No batch
#' correction happens here - the resulting object's own (un-integrated)
#' PCA/UMAP downstream is exactly the baseline an integration method gets
#' compared against.
merge_for_integration <- function(object_list) {
  objs <- unname(object_list)
  merged <- if (length(objs) == 2) merge(objs[[1]], y = objs[[2]]) else merge(objs[[1]], y = objs[-1])
  stopifnot(!anyDuplicated(colnames(merged)))
  merged
}

#' Adds a harmonized `major_population` metadata column using
#' metadata/label_mapping.tsv, which maps A's broad_label alone, and B's
#' (ano.l1, ano.l2) PAIR (B has documented cross_population_conflict cases
#' where the same ano.l1 maps to different populations depending on ano.l2 -
#' see DATA_DICTIONARY.md), to a shared population category. Cells with no
#' match (or mapped to the mapping's own "Unresolved" bucket) get
#' major_population = "Unresolved" rather than being silently dropped, since
#' where those cells land is itself informative for a qualitative check.
add_major_population <- function(obj, label_mapping) {
  meta <- obj@meta.data
  map_a <- label_mapping[label_mapping$cohort == "A", c("broad_label", "population")]
  map_b <- label_mapping[label_mapping$cohort == "B", c("broad_label", "detailed_label", "population")]

  pop <- rep(NA_character_, nrow(meta))
  is_a <- meta$cohort == "A"
  is_b <- meta$cohort == "B"

  pop[is_a] <- map_a$population[match(meta$broad_label[is_a], map_a$broad_label)]
  b_key <- paste(meta$ano.l1[is_b], meta$ano.l2[is_b], sep = "\r")
  map_b_key <- paste(map_b$broad_label, map_b$detailed_label, sep = "\r")
  pop[is_b] <- map_b$population[match(b_key, map_b_key)]

  pop[is.na(pop)] <- "Unresolved"
  obj$major_population <- pop
  n_unresolved <- sum(pop == "Unresolved")
  cat(sprintf("  major_population: %d/%d cells (%.1f%%) Unresolved/unmatched\n",
              n_unresolved, length(pop), 100 * n_unresolved / length(pop)))
  obj
}
