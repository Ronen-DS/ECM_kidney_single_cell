## ECM panel expression across major cell populations, split into 2
## side-by-side panels (Tumor, Adjacent) - the second robustness axis
## (tumor versus adjacent tissue), pooling across all three cohorts within
## each panel (cohort-level splitting was done separately in
## ecm_panel_by_population_per_cohort.R).
##
## Relative expression (z) is computed ONCE per gene across ALL
## population x condition groups together (not separately within each
## condition) - unlike the cohort-split plot, where within-cohort scaling is
## the right choice because absolute expression levels genuinely differ
## between nuclei (A), whole cells (B), and C's derivative counts (a
## technical/platform difference). Tumor and Adjacent share the same
## measurement scale within each cohort, so z-scoring them separately would
## make the two panels' colors incomparable - exactly the comparison this
## plot exists to show. A population x condition group with fewer than 20
## cells is marked with a black outline rather than hidden, since the data
## point itself is still worth seeing, just flagged as unreliable.
##
## Input:  merged_ABC_annotated.rds, seurat_{A,B,C}_normalized.rds
## Output: results/tables/ecm_panel_by_population_per_condition.tsv
##         results/figures/ecm_panel_by_population_per_condition.png
##
## Run: Rscript R/ecm_panel_by_population_per_condition.R

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
library(Matrix)
library(ggplot2)

ecm_genes <- c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1",
               "COL6A2", "FN1", "SPARC", "DCN", "LUM")

cat("Loading merged_ABC_annotated.rds for metadata only...\n")
obj_meta <- readRDS(file.path(SCRATCH_DIR, "merged_ABC_annotated.rds"))
meta <- data.frame(
  cell_key = colnames(obj_meta),
  cohort = obj_meta$cohort,
  condition = obj_meta$condition,
  cell_annotation = obj_meta$cell_annotation,
  stringsAsFactors = FALSE
)
rm(obj_meta)
gc(verbose = FALSE)
cat(sprintf("  %d cells\n", nrow(meta)))
cat("Condition counts:\n")
print(table(meta$condition))

ecm_frames <- list()
for (cohort in c("A", "B", "C")) {
  cat(sprintf("Cohort %s...\n", cohort))
  obj <- readRDS(file.path(SCRATCH_DIR, sprintf("seurat_%s_normalized.rds", cohort)))
  mat <- GetAssayData(obj, assay = "RNA", layer = "data")[ecm_genes, , drop = FALSE]
  df <- as.data.frame(as.matrix(t(mat)))
  df$cell_key <- colnames(obj)
  ecm_frames[[cohort]] <- df[, c("cell_key", ecm_genes)]
  rm(obj, mat)
  gc(verbose = FALSE)
}
ecm_df <- do.call(rbind, ecm_frames)
combined <- merge(meta, ecm_df, by = "cell_key")
cat(sprintf("Combined: %d cells\n", nrow(combined)))

## --- Per-(condition, population) summary, raw (no z-score yet) --------------
MIN_CELLS <- 20

summarize_group <- function(data, ecm_genes, cond) {
  populations <- sort(unique(data$cell_annotation))
  rows <- list()
  for (p in populations) {
    sub <- data[data$cell_annotation == p, ]
    n_cells <- nrow(sub)
    if (n_cells == 0) next
    for (g in ecm_genes) {
      vals <- sub[[g]]
      mean_expr <- mean(vals)
      pct_expr <- 100 * sum(vals > 0) / n_cells
      rows[[length(rows) + 1]] <- data.frame(condition = cond, population = p, gene = g,
                                              mean_expr = mean_expr, pct_expr = pct_expr, n_cells = n_cells)
    }
  }
  do.call(rbind, rows)
}

all_rows <- list()
for (cond in c("Tumor", "Adjacent")) {
  sub <- combined[combined$condition == cond, ]
  tbl <- summarize_group(sub, ecm_genes, cond)
  all_rows[[cond]] <- tbl
  cat(sprintf("\n%s: %d cells, %d populations with cells present\n", cond, nrow(sub), length(unique(tbl$population))))
  missing_pop <- setdiff(sort(unique(combined$cell_annotation)), unique(tbl$population))
  if (length(missing_pop) > 0) cat(sprintf("  No cells in: %s\n", paste(missing_pop, collapse = ", ")))
}
full_tbl <- do.call(rbind, all_rows)

## Relative expression: z-score each gene ONCE across all (condition,
## population) groups together, so colors are comparable between the two
## panels - see header comment.
full_tbl$expr_z <- ave(full_tbl$mean_expr, full_tbl$gene, FUN = function(x) {
  if (length(x) < 2 || sd(x) == 0) return(rep(NA_real_, length(x)))
  (x - mean(x)) / sd(x)
})
full_tbl$low_n <- full_tbl$n_cells < MIN_CELLS
cat(sprintf("\n%d/%d (condition, population) groups have < %d cells (marked, not dropped)\n",
            sum(full_tbl$low_n) / length(ecm_genes), nrow(full_tbl) / length(ecm_genes), MIN_CELLS))

tables_dir <- file.path(repo, "results", "tables")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
write.table(full_tbl, file.path(tables_dir, "ecm_panel_by_population_per_condition.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

## --- 2-panel figure ----------------------------------------------------------
fig_dir <- file.path(repo, "results", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

full_tbl$gene <- factor(full_tbl$gene, levels = rev(ecm_genes))
full_tbl$condition <- factor(full_tbl$condition, levels = c("Tumor", "Adjacent"))
full_tbl$low_n_label <- ifelse(full_tbl$low_n, "< 20 cells", ">= 20 cells")
p <- ggplot(full_tbl[!is.na(full_tbl$mean_expr), ],
            aes(x = population, y = gene, size = pct_expr, fill = expr_z, color = low_n_label)) +
  geom_point(shape = 21, stroke = 1.1) +
  facet_wrap(~condition, nrow = 1) +
  scale_fill_gradient2(low = "blue", mid = "grey90", high = "red", midpoint = 0, name = "Relative\nexpression (z)") +
  scale_color_manual(values = c(">= 20 cells" = "white", "< 20 cells" = "black"), name = "Group size",
                      guide = guide_legend(override.aes = list(fill = "grey70", size = 5))) +
  scale_size_continuous(name = "% expressing", range = c(0, 10)) +
  labs(title = "ECM panel expression by population, split by condition (Tumor vs Adjacent)",
       subtitle = "Relative expression z-scored once across both panels together - colors ARE comparable between Tumor and Adjacent",
       x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        strip.text = element_text(face = "bold", size = 12))
ggsave(file.path(fig_dir, "ecm_panel_by_population_per_condition.png"), p, width = 12, height = 6.5, dpi = 150, bg = "white")

cat("\nWrote results/figures/ecm_panel_by_population_per_condition.png\n")
cat("Done.\n")
