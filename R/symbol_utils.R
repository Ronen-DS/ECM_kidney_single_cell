## Helpers for collapsing Ensembl-ID duplicate rows down to one row per gene
## SYMBOL, by summing counts. Needed because cohort B's matrix is already
## symbol-indexed (no Ensembl IDs at all) - to eventually put A/B/C on a
## shared gene-identifier space, A and C's Ensembl-ID rows need to become
## symbol rows too, and a small number of symbols in A/C map to >1 distinct
## Ensembl ID (GENCODE reference quirks - see the duplicate-ID audit) so
## those need summing rather than an arbitrary pick.

library(Matrix)

#' Sums duplicate-symbol rows of a genes x cells count matrix into one row
#' per unique symbol. `features_df` must have one row per row of `counts`,
#' in the same order, with an id column matching rownames(counts) and a
#' symbol column giving the target grouping. Implemented as a sparse
#' indicator-matrix product (crossprod) rather than base::rowsum(), so it
#' stays sparse and scales to full-size count matrices.
collapse_symbol_duplicates <- function(counts, features_df,
                                        id_col = "gene_id", symbol_col = "gene_symbol") {
  stopifnot(identical(rownames(counts), features_df[[id_col]]))
  symbols <- features_df[[symbol_col]]
  unique_symbols <- unique(symbols)
  indicator <- Matrix::sparseMatrix(
    i = seq_along(symbols),
    j = match(symbols, unique_symbols),
    x = 1,
    dims = c(length(symbols), length(unique_symbols))
  )
  new_counts <- Matrix::crossprod(indicator, counts)
  rownames(new_counts) <- unique_symbols
  methods::as(new_counts, "CsparseMatrix")
}

#' A (cohort, gene_symbol, n_duplicates, ensembl_ids) row for every symbol
#' mapped to more than one Ensembl ID in `features_df` - the set actually
#' summed together by collapse_symbol_duplicates() for that cohort.
duplicate_symbol_table <- function(features_df, cohort_label,
                                    id_col = "gene_id", symbol_col = "gene_symbol") {
  by_symbol <- split(features_df[[id_col]], features_df[[symbol_col]])
  dup <- by_symbol[lengths(by_symbol) > 1]
  if (length(dup) == 0) {
    return(data.frame(cohort = character(0), gene_symbol = character(0),
                       n_duplicates = integer(0), ensembl_ids = character(0),
                       stringsAsFactors = FALSE))
  }
  data.frame(
    cohort = cohort_label,
    gene_symbol = names(dup),
    n_duplicates = lengths(dup),
    ensembl_ids = vapply(dup, paste, character(1), collapse = "; "),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
