## Central configuration - paths, R library, and R version pin.
##
## DATA_ROOT is intentionally NOT hardcoded here (the applicant-package data
## folder is not part of this git repo - "reference the input manifest
## instead of committing the data"). Set it once per machine via an
## environment variable, e.g. in a local (gitignored) .Renviron file at the
## repo root:
##
##   CCRCC_DATA_ROOT=C:/path/to/applicant-package
##
## or by calling Sys.setenv(CCRCC_DATA_ROOT = "...") before sourcing this
## file. See README.md for the value used to produce the reported results.

## Package library location(s). Packages for this project were installed into
## a per-user library (not the default R installation library, which needs
## admin rights on Windows) - see README.md "Environment" section for exactly
## which R build/library this was run against. Adding this path is a no-op on
## any machine where it doesn't exist.
user_lib <- file.path(Sys.getenv("USERPROFILE", unset = Sys.getenv("HOME")), "R", "library-x64")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

DATA_ROOT <- Sys.getenv("CCRCC_DATA_ROOT", unset = NA)
if (is.na(DATA_ROOT) || !nzchar(DATA_ROOT)) {
  stop(
    "CCRCC_DATA_ROOT is not set. Point it at the applicant-package folder, ",
    "e.g. Sys.setenv(CCRCC_DATA_ROOT = '/path/to/applicant-package') or via ",
    "a local .Renviron entry - see README.md."
  )
}
DATA_ROOT <- normalizePath(DATA_ROOT, mustWork = TRUE)

RNA_DIR <- file.path(DATA_ROOT, "data", "rna")
BULK_DIR <- file.path(DATA_ROOT, "data", "bulk")
METADATA_DIR <- file.path(DATA_ROOT, "metadata")

COHORT_A_DIR <- file.path(RNA_DIR, "cohort_A_GSE227898")
COHORT_B_DIR <- file.path(RNA_DIR, "cohort_B_GSE178481")
COHORT_C_DIR <- file.path(RNA_DIR, "cohort_C_GSE159115")

SAMPLE_MANIFEST_PATH <- file.path(METADATA_DIR, "sample_manifest.tsv")
ECM_PANEL_PATH <- file.path(METADATA_DIR, "ecm_program.tsv")
LABEL_MAPPING_PATH <- file.path(METADATA_DIR, "label_mapping.tsv")
MULTIMODAL_CROSSWALK_PATH <- file.path(METADATA_DIR, "multimodal_crosswalk.tsv")
PAPER_POPULATION_COUNTS_PATH <- file.path(METADATA_DIR, "paper_population_counts.tsv")
BULK_SAMPLES_PATH <- file.path(METADATA_DIR, "bulk_samples.tsv")
BULK_COLUMN_MANIFEST_PATH <- file.path(METADATA_DIR, "bulk_column_manifest.tsv")

ANCHOR_ANNOTATIONS_PATH <- file.path(METADATA_DIR, "original", "anchor_annotations.tsv.gz")
ANCHOR_SAMPLES_PATH <- file.path(METADATA_DIR, "original", "anchor_samples.tsv")
COHORT_B_ANNOTATIONS_PATH <- file.path(METADATA_DIR, "original", "cohort_b_annotations.csv")
COHORT_C_SAMPLES_PATH <- file.path(METADATA_DIR, "original", "cohort_C_samples.tsv")

## A persistent working/cache directory for extracted archives and expensive
## intermediate results (e.g. the full cohort A load, ~17 min to rebuild) -
## kept outside the repo and outside DATA_ROOT, so the read-only input
## package is never modified. Deliberately NOT tempdir(): that's unique to
## each R session, so anything cached there would be silently orphaned and
## invisible to the next script run. tools::R_user_dir() is the standard,
## portable per-user cache location (works the same way on any machine).
SCRATCH_DIR <- tools::R_user_dir("ccrcc_ecm_analysis", "cache")
dir.create(SCRATCH_DIR, showWarnings = FALSE, recursive = TRUE)

RANDOM_SEED <- 1L
