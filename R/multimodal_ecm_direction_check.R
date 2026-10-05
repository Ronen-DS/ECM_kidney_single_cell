## Part 3 (continued): compares the Tumor-vs-Adjacent DIRECTION of the chosen
## ECM finding (Endothelial SPARC/FN1/panel-level up-in-Tumor signal,
## established in R/ecm_endothelial_models.R - strongest/most significant in
## cohort B, null in cohort C, cohort A alone was never significant on its
## own) against linked bulk RNA and proteomics evidence, within cohort A only
## (the only cohort with bulk/protein links - see multimodal_matched_subset.R).
##
## This is a DIRECTION/plausibility check, not a validation of cell-type
## specificity: bulk RNA and proteomics are WHOLE-TISSUE measurements (all
## cell types mixed), so a concordant bulk signal is also consistent with a
## shift in overall tissue composition (e.g. more stromal content in Tumor),
## not proof the effect originates in Endothelial cells specifically. No
## specimen-level RNA-protein correlation or joint embedding is computed here
## (explicitly out of scope per the task).
##
## Per-assay valid sample sets differ (see multimodal_matched_subset.R):
##   bulk RNA:  all crosswalk rows with a bulk_sample_id (matching_status !=
##              "no_author_bulk_link") - 29 specimens, since a bulk RNA link
##              only requires donor+condition matching, not a protein
##              aliquot's parent-specimen status.
##   proteomics: only "same_parent_specimen" rows - 27 specimens.
## Donor pairing (both Tumor AND Adjacent valid) is therefore assay-specific
## and computed separately for each.
##
## Panel genes are log2, already-normalized values in both bulk files, so a
## "panel score" is the MEAN of the 10 genes' log2 values per sample (NOT a
## sum of raw counts - that pseudobulk-specific method doesn't apply to
## already-log2 bulk/protein data). Any panel gene absent from a file's gene
## space is excluded from that file's mean for ALL samples and reported, not
## zero-filled.
##
## Input:  metadata/multimodal_crosswalk.tsv,
##         data/bulk/CCRCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_{Normal,Tumor}.txt,
##         data/bulk/CCRCC_proteomics_gene_abundance_log2_reference_intensity_normalized_{Normal,Tumor}.txt,
##         <SCRATCH_DIR>/seurat_A_raw.rds (misc$features: gene_id<->gene_symbol lookup)
## Output: results/multimodal/ecm_direction_bulk_vs_protein.tsv (long table)
##         results/multimodal/ecm_direction_paired_summary.tsv
##
## Run: Rscript R/multimodal_ecm_direction_check.R

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

out_dir <- file.path(repo, "results", "multimodal")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ecm_genes <- c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1",
               "COL6A2", "FN1", "SPARC", "DCN", "LUM")

cat("Loading cohort A gene_id<->gene_symbol lookup (seurat_A_raw.rds@misc$features)...\n")
obj <- readRDS(file.path(SCRATCH_DIR, "seurat_A_raw.rds"))
feat <- obj@misc$features
rm(obj); gc(verbose = FALSE)
id_map <- setNames(feat$gene_id, feat$gene_symbol)
stopifnot(all(ecm_genes %in% names(id_map)))  # confirmed earlier: all 10 present in cohort A's reference

cw <- read.delim(MULTIMODAL_CROSSWALK_PATH, stringsAsFactors = FALSE)

bulk_valid <- unique(cw[cw$matching_status != "no_author_bulk_link", c("donor_id", "condition")])
protein_valid <- unique(cw[cw$matching_status == "same_parent_specimen", c("donor_id", "condition")])
cat(sprintf("Bulk RNA valid specimens: %d (donor x condition)\n", nrow(bulk_valid)))
cat(sprintf("Protein valid specimens: %d (donor x condition)\n", nrow(protein_valid)))

extract_panel <- function(file_path, valid_df, ecm_genes, id_map) {
  cat(sprintf("  Reading %s...\n", basename(file_path)))
  mat <- read.delim(file_path, row.names = 1, check.names = FALSE)
  ## Stripping the Ensembl version suffix creates a handful of duplicate base
  ## IDs elsewhere in the genome (PAR-region/retired-ID duplicates, confirmed
  ## none of the 10 ECM panel genes are among them) - match positionally via
  ## match() (first hit) instead of reassigning rownames globally, which would
  ## require uniqueness across the WHOLE file.
  base_ids <- sub("\\..*$", "", rownames(mat))

  present_genes <- ecm_genes[id_map[ecm_genes] %in% base_ids]
  missing_genes <- setdiff(ecm_genes, present_genes)
  if (length(missing_genes) > 0) cat(sprintf("    Panel genes ABSENT from this file (excluded from panel mean, not zero-filled): %s\n", paste(missing_genes, collapse = ", ")))

  donors_in_file <- intersect(valid_df$donor_id, colnames(mat))
  missing_donors <- setdiff(valid_df$donor_id, colnames(mat))
  if (length(missing_donors) > 0) cat(sprintf("    NOTE: donor(s) in crosswalk but not a column in this file: %s\n", paste(missing_donors, collapse = ", ")))

  row_idx <- match(id_map[present_genes], base_ids)
  stopifnot(!anyNA(row_idx))
  sub_mat <- mat[row_idx, donors_in_file, drop = FALSE]
  rownames(sub_mat) <- present_genes
  list(values = sub_mat, present_genes = present_genes, missing_genes = missing_genes)
}

build_assay_table <- function(assay_label, tumor_file, normal_file, valid_df) {
  cat(sprintf("\n=== %s ===\n", assay_label))
  tumor_rows <- valid_df[valid_df$condition == "Tumor", ]
  adj_rows <- valid_df[valid_df$condition == "Adjacent", ]

  tumor_ext <- extract_panel(tumor_file, tumor_rows, ecm_genes, id_map)
  adj_ext <- extract_panel(normal_file, adj_rows, ecm_genes, id_map)
  present_genes <- intersect(tumor_ext$present_genes, adj_ext$present_genes)
  cat(sprintf("  Panel genes used in both Tumor and Adjacent files: %s (n=%d)\n",
              paste(present_genes, collapse = ", "), length(present_genes)))

  rows <- list()
  for (donor in tumor_rows$donor_id[tumor_rows$donor_id %in% colnames(tumor_ext$values)]) {
    v <- tumor_ext$values[present_genes, donor]
    rows[[length(rows) + 1]] <- data.frame(assay = assay_label, donor_id = donor, condition = "Tumor",
                                            SPARC = tumor_ext$values["SPARC", donor], FN1 = tumor_ext$values["FN1", donor],
                                            panel_mean = mean(v))
  }
  for (donor in adj_rows$donor_id[adj_rows$donor_id %in% colnames(adj_ext$values)]) {
    v <- adj_ext$values[present_genes, donor]
    rows[[length(rows) + 1]] <- data.frame(assay = assay_label, donor_id = donor, condition = "Adjacent",
                                            SPARC = adj_ext$values["SPARC", donor], FN1 = adj_ext$values["FN1", donor],
                                            panel_mean = mean(v))
  }
  do.call(rbind, rows)
}

bulk_tbl <- build_assay_table("bulk_RNA",
                               file.path(DATA_ROOT, "data", "bulk", "CCRCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt"),
                               file.path(DATA_ROOT, "data", "bulk", "CCRCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt"),
                               bulk_valid)
protein_tbl <- build_assay_table("protein",
                                  file.path(DATA_ROOT, "data", "bulk", "CCRCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"),
                                  file.path(DATA_ROOT, "data", "bulk", "CCRCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt"),
                                  protein_valid)

all_tbl <- rbind(bulk_tbl, protein_tbl)
write.table(all_tbl, file.path(out_dir, "ecm_direction_bulk_vs_protein.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## --- Endothelial snRNA pseudobulk coverage lookup, for cross-referencing ---
## whether a bulk/protein-paired donor's Tumor-vs-Adjacent whole-tissue
## difference is backed by its OWN paired Endothelial single-cell data (the
## cell type the original ECM finding is about), or is generic whole-tissue
## evidence with no single-cell Endothelial comparison to cross-check against
## in that same donor (one or both conditions never reached the >=20-cell
## Endothelial pseudobulk threshold - see ecm_endothelial_pseudobulk_build.R).
endo_pb <- readRDS(file.path(SCRATCH_DIR, "pseudobulk_endothelial.rds"))$sample_table
endo_pb_a <- endo_pb[endo_pb$cohort == "A", ]
endothelial_coverage <- function(donor) {
  sub <- endo_pb_a[endo_pb_a$donor_id == donor, ]
  n_tumor <- sub$n_cells[sub$condition == "Tumor"]
  n_adjacent <- sub$n_cells[sub$condition == "Adjacent"]
  data.frame(
    endo_n_cells_tumor = if (length(n_tumor) == 1) n_tumor else 0L,
    endo_n_cells_adjacent = if (length(n_adjacent) == 1) n_adjacent else 0L,
    has_paired_endothelial_snRNA = length(n_tumor) == 1 && length(n_adjacent) == 1
  )
}

## --- Paired donors per assay: Tumor - Adjacent difference ------------------
summarize_assay <- function(tbl, assay_label) {
  sub <- tbl[tbl$assay == assay_label, ]
  key <- table(sub$donor_id)
  paired_donors <- names(key)[key == 2]
  unpaired_donors <- names(key)[key == 1]
  cat(sprintf("\n%s: %d donors total (%d paired, %d unpaired)\n", assay_label,
              length(key), length(paired_donors), length(unpaired_donors)))

  diff_rows <- list()
  for (donor in paired_donors) {
    t_row <- sub[sub$donor_id == donor & sub$condition == "Tumor", ]
    a_row <- sub[sub$donor_id == donor & sub$condition == "Adjacent", ]
    diff_rows[[donor]] <- cbind(
      data.frame(assay = assay_label, donor_id = donor,
                 SPARC_diff = t_row$SPARC - a_row$SPARC,
                 FN1_diff = t_row$FN1 - a_row$FN1,
                 panel_mean_diff = t_row$panel_mean - a_row$panel_mean),
      endothelial_coverage(donor)
    )
  }
  diffs <- if (length(diff_rows) > 0) do.call(rbind, diff_rows) else
    data.frame(assay = character(0), donor_id = character(0), SPARC_diff = numeric(0), FN1_diff = numeric(0), panel_mean_diff = numeric(0),
               endo_n_cells_tumor = integer(0), endo_n_cells_adjacent = integer(0), has_paired_endothelial_snRNA = logical(0))

  cat("  Paired donor differences (Tumor - Adjacent), with Endothelial snRNA coverage:\n")
  print(diffs, row.names = FALSE)
  cat(sprintf("  %d/%d paired donors have their OWN paired Endothelial snRNA data (>=20 cells in both conditions) to cross-reference against.\n",
              sum(diffs$has_paired_endothelial_snRNA), nrow(diffs)))

  if (nrow(diffs) >= 2) {
    for (feat_name in c("SPARC_diff", "FN1_diff", "panel_mean_diff")) {
      tt <- t.test(diffs[[feat_name]])
      cat(sprintf("  Paired t-test %s: mean diff=%.3f, 95%% CI (%.3f, %.3f), P=%.3f (n=%d - exploratory only, not a formal inference given sample size)\n",
                  feat_name, tt$estimate, tt$conf.int[1], tt$conf.int[2], tt$p.value, nrow(diffs)))
    }
  }
  cat("  Unpaired donors (single condition only, Tumor-vs-Adjacent group means shown for context):\n")
  unpaired_sub <- sub[sub$donor_id %in% unpaired_donors, ]
  print(aggregate(cbind(SPARC, FN1, panel_mean) ~ condition, data = unpaired_sub, FUN = mean), row.names = FALSE)

  diffs
}

bulk_diffs <- summarize_assay(all_tbl, "bulk_RNA")
protein_diffs <- summarize_assay(all_tbl, "protein")

paired_summary <- rbind(bulk_diffs, protein_diffs)
write.table(paired_summary, file.path(out_dir, "ecm_direction_paired_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== Reference: Endothelial snRNA pseudobulk direction (from ecm_endothelial_models.R) ===\n")
cat("  A only:  panel-sum logFC +0.11 (P=0.83, NOT significant); SPARC/FN1 not individually significant\n")
cat("  B only:  panel-sum logFC +1.62 (P=9.1e-8, significant); SPARC +1.89, FN1 +1.09 (both significant)\n")
cat("  C only:  panel-sum logFC +0.03 (P=0.97, NOT significant, null)\n")
cat("  -> The 'chosen ECM finding' (Tumor-up SPARC/FN1/panel) is a cohort-B-driven\n")
cat("     signal; cohort A's OWN snRNA data never showed it significantly. Bulk RNA\n")
cat("     and protein data exist ONLY for cohort A donors, so this comparison tests\n")
cat("     generalization to a cohort/compartment where the single-cell signal was\n")
cat("     already null - temper conclusions accordingly.\n")

cat("\nDone.\n")
