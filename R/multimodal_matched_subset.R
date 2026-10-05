## Part 3: validates donor/specimen/aliquot correspondence across modalities
## using the supplied crosswalk (metadata/multimodal_crosswalk.tsv), and
## reports the matched subset usable for a direct snRNA-bulkRNA-protein
## comparison, with exclusions/ambiguity explained.
##
## Scope: the crosswalk only covers cohort A. Cohorts B and C use entirely
## separate donor-ID namespaces (RCC-PRx / "Patient SS_xxxx") from distinct
## GEO studies with no bulk RNA or proteomics companion data in this package -
## confirmed below by checking the crosswalk's cohort column directly, rather
## than assumed. The multi-omic comparison in Part 3 is therefore necessarily
## cohort-A-only.
##
## Link quality is NOT uniform across modalities for a given row:
##  - snRNA <-> bulk RNA: matched by donor_id + condition column only
##    ("bulk_RNA_provenance" = "donor_condition_column; exact_RNA_aliquot_not_
##    reconciled" on EVERY row that has a bulk_sample_id) - i.e. bulk RNA is
##    assumed to be the one bulk sample for that donor+condition, but was
##    never verified at the aliquot/parent-specimen level. This caveat
##    applies even to "same_parent_specimen" rows.
##  - snRNA <-> proteomics: matched by shared PARENT SPECIMEN ID (aliquot-
##    level, author-reconciled) - this is what `matching_status` /
##    `primary_specimen_link` actually certify.
##
## matching_status categories found:
##  - same_parent_specimen (primary_specimen_link=True): clean 3-way match -
##    snRNA aliquot's parent specimen == protein aliquot's parent specimen,
##    AND a donor+condition-matched bulk RNA sample exists.
##  - no_author_bulk_link: neither bulk RNA nor protein data exists for this
##    SPECIFIC snRNA aliquot (though another aliquot from the same donor+
##    condition may be matched instead - checked explicitly below).
##  - different_parent_specimen: the protein aliquot comes from a physically
##    different parent specimen than the snRNA aliquot - excluded from any
##    snRNA<->protein comparison (bulk RNA, being only donor+condition
##    matched in the first place, is not inherently invalidated by this).
##  - multiple_protein_specimen_ids: the protein data maps ambiguously to
##    >1 parent specimen - which one corresponds to this snRNA aliquot can't
##    be resolved, so excluded from snRNA<->protein comparison.
##
## Output: results/multimodal/matched_subset.tsv (clean 3-way matches)
##         results/multimodal/excluded_specimens.tsv (with reasons)
##         results/multimodal/donor_condition_coverage.tsv (per donor x
##           condition: how many aliquots, how many matched, net coverage)
##
## Run: Rscript R/multimodal_matched_subset.R

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

cw <- read.delim(MULTIMODAL_CROSSWALK_PATH, stringsAsFactors = FALSE)
cat(sprintf("Crosswalk: %d rows\n", nrow(cw)))
cat("Cohorts present in crosswalk: "); print(table(cw$cohort))
stopifnot(all(cw$cohort == "A"))  # confirms B/C have no bulk/proteomics linkage in this package

cat("\n=== matching_status breakdown ===\n")
print(table(cw$matching_status))
cat("\nprimary_specimen_link breakdown:\n")
print(table(cw$primary_specimen_link))

n_donors <- length(unique(cw$donor_id))
n_donor_conditions <- length(unique(paste(cw$donor_id, cw$condition)))
cat(sprintf("\nCohort A: %d snRNA specimens, %d unique donors, %d unique donor x condition combinations\n",
            nrow(cw), n_donors, n_donor_conditions))

## --- Matched subset: clean 3-way link -----------------------------------
matched <- cw[cw$matching_status == "same_parent_specimen", ]
write.table(matched, file.path(out_dir, "matched_subset.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("\nMatched subset (snRNA + bulk RNA + protein, same parent specimen): %d/%d specimens (%.1f%%)\n",
            nrow(matched), nrow(cw), 100 * nrow(matched) / nrow(cw)))
cat(sprintf("  Spanning %d unique donors, %d donor x condition combinations\n",
            length(unique(matched$donor_id)), length(unique(paste(matched$donor_id, matched$condition)))))

## --- Excluded specimens, with explicit reasons --------------------------
excluded <- cw[cw$matching_status != "same_parent_specimen", ]
excluded$exclusion_reason <- ifelse(
  excluded$matching_status == "no_author_bulk_link",
  "No bulk RNA or protein data for this specific snRNA aliquot (author never profiled it at the bulk/protein level)",
  ifelse(excluded$matching_status == "different_parent_specimen",
         "Protein aliquot is from a DIFFERENT parent specimen than this snRNA aliquot - not the same physical tissue piece; exclude from snRNA<->protein comparison specifically",
         ifelse(excluded$matching_status == "multiple_protein_specimen_ids",
                "Protein data maps ambiguously to >1 parent specimen - cannot determine which corresponds to this snRNA aliquot; exclude from snRNA<->protein comparison",
                "unclassified")))

## For no_author_bulk_link rows: check whether the SAME donor+condition has
## a different, matched aliquot covering it (a redundant-replicate case) or
## whether this is the donor+condition's only aliquot (a full coverage gap).
excluded$donor_condition <- paste(excluded$donor_id, excluded$condition)
matched_donor_conditions <- unique(paste(matched$donor_id, matched$condition))
excluded$covered_by_other_aliquot <- excluded$donor_condition %in% matched_donor_conditions

write.table(excluded[, c("sample_id", "donor_id", "condition", "matching_status", "exclusion_reason", "covered_by_other_aliquot")],
            file.path(out_dir, "excluded_specimens.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("\nExcluded specimens: %d/%d\n", nrow(excluded), nrow(cw)))
print(excluded[, c("sample_id", "donor_id", "condition", "matching_status", "covered_by_other_aliquot")], row.names = FALSE)

cat("\n--- Exclusion reason breakdown ---\n")
print(table(excluded$matching_status))
cat(sprintf("\nOf the excluded specimens, %d/%d belong to a donor x condition that IS otherwise covered by a different matched aliquot (redundant replicate, no net coverage loss).\n",
            sum(excluded$covered_by_other_aliquot), nrow(excluded)))
cat(sprintf("The remaining %d exclusion(s) are the ONLY aliquot for their donor x condition - that donor x condition has NO usable snRNA<->protein link at all:\n",
            sum(!excluded$covered_by_other_aliquot)))
print(excluded[!excluded$covered_by_other_aliquot, c("sample_id", "donor_id", "condition", "matching_status")], row.names = FALSE)

## --- Donor x condition coverage summary ---------------------------------
cov <- aggregate(sample_id ~ donor_id + condition, data = cw, FUN = length)
names(cov)[3] <- "n_aliquots"
cov_matched <- aggregate(sample_id ~ donor_id + condition, data = matched, FUN = length)
names(cov_matched)[3] <- "n_matched"
cov <- merge(cov, cov_matched, by = c("donor_id", "condition"), all.x = TRUE)
cov$n_matched[is.na(cov$n_matched)] <- 0
cov$fully_unmatched <- cov$n_matched == 0
cov <- cov[order(cov$donor_id, cov$condition), ]
write.table(cov, file.path(out_dir, "donor_condition_coverage.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("\n=== Donor x condition coverage: %d combinations total, %d with zero matched aliquot ===\n",
            nrow(cov), sum(cov$fully_unmatched)))
if (sum(cov$fully_unmatched) > 0) print(cov[cov$fully_unmatched, ], row.names = FALSE)

cat("\nDone.\n")
