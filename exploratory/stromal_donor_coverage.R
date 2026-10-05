## Coverage table for the donor-level robustness axis: how many
## Stromal_fibro_myo cells does each donor x condition x cohort cell have,
## and how many donors clear plausible minimum-cell thresholds (10, 20) in
## BOTH Tumor and Adjacent simultaneously (required for a paired/matched
## donor-level comparison).
##
## Input:  merged_ABC_annotated.rds
## Output: results/tables/stromal_donor_coverage.tsv (full donor x condition x cohort table)
##
## Run: Rscript R/stromal_donor_coverage.R

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

cat("Loading merged_ABC_annotated.rds for metadata only...\n")
obj <- readRDS(file.path(SCRATCH_DIR, "merged_ABC_annotated.rds"))
meta <- obj@meta.data
rm(obj); gc(verbose = FALSE)

stopifnot(all(c("donor_id", "condition", "cohort", "cell_annotation") %in% colnames(meta)))

stromal <- meta[meta$cell_annotation == "Stromal_fibro_myo", ]
cat(sprintf("Total Stromal_fibro_myo cells: %d\n", nrow(stromal)))
cat("By cohort:\n"); print(table(stromal$cohort))

## --- Full donor x condition x cohort table (all donors in the data, 0-filled) ---
all_donors <- unique(meta[, c("donor_id", "cohort", "condition")])
counts <- aggregate(rep(1, nrow(stromal)),
                     by = list(donor_id = stromal$donor_id,
                               cohort = stromal$cohort,
                               condition = stromal$condition),
                     FUN = sum)
names(counts)[4] <- "n_cells"

full_tbl <- merge(all_donors, counts, by = c("donor_id", "cohort", "condition"), all.x = TRUE)
full_tbl$n_cells[is.na(full_tbl$n_cells)] <- 0
full_tbl <- full_tbl[order(full_tbl$cohort, full_tbl$donor_id, full_tbl$condition), ]

tables_dir <- file.path(repo, "results", "tables")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
write.table(full_tbl, file.path(tables_dir, "stromal_donor_coverage.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("\nWrote %s (%d donor x condition x cohort rows)\n",
            file.path(tables_dir, "stromal_donor_coverage.tsv"), nrow(full_tbl)))

## --- Wide view: Tumor vs Adjacent cell counts per donor, for donors that have BOTH conditions present at all ---
wide <- reshape(full_tbl, idvar = c("donor_id", "cohort"), timevar = "condition", direction = "wide")
names(wide) <- sub("^n_cells\\.", "", names(wide))
if (!"Tumor" %in% names(wide)) wide$Tumor <- 0
if (!"Adjacent" %in% names(wide)) wide$Adjacent <- 0
wide$Tumor[is.na(wide$Tumor)] <- 0
wide$Adjacent[is.na(wide$Adjacent)] <- 0
wide <- wide[order(wide$cohort, wide$donor_id), c("donor_id", "cohort", "Tumor", "Adjacent")]

cat("\nPer-donor Stromal_fibro_myo cell counts (Tumor vs Adjacent):\n")
print(wide, row.names = FALSE)

## --- Threshold check: donors with >= N cells in BOTH Tumor AND Adjacent ---
for (thresh in c(10, 20)) {
  qualifying <- wide[wide$Tumor >= thresh & wide$Adjacent >= thresh, ]
  cat(sprintf("\nDonors with >= %d Stromal_fibro_myo cells in BOTH Tumor and Adjacent: %d / %d donors with any Stromal_fibro_myo cells\n",
              thresh, nrow(qualifying), sum(wide$Tumor > 0 | wide$Adjacent > 0)))
  if (nrow(qualifying) > 0) {
    cat("  Qualifying donors:\n")
    print(qualifying, row.names = FALSE)
  }
  cat("  By cohort:\n")
  print(table(qualifying$cohort))
}

cat("\nDone.\n")
