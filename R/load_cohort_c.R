## Cohort C (GSE159115): legacy 10x HDF5, one file per specimen, under a
## GRCh38 group holding a CSC sparse matrix (data/indices/indptr/shape) plus
## genes (Ensembl ID) and gene_names (symbol). No cell labels are supplied -
## every barcode in the file is already the deposited/retained cell set (no
## additional filtering or join needed, unlike cohort A's raw droplets).

library(hdf5r)
library(Matrix)

#' Reads one cohort C .h5 file as a sparse dgCMatrix (genes x cells), plus a
#' features data.frame with gene_id/gene_symbol (mirrors cohort A's
#' features table shape for later cross-cohort harmonization).
load_cohort_c_sample <- function(sample_id, h5_path) {
  f <- H5File$new(h5_path, mode = "r")
  on.exit(f$close_all())
  g <- f[["GRCh38"]]

  shape <- g[["shape"]]$read() # c(n_genes, n_cells)
  mat <- Matrix::sparseMatrix(
    i = g[["indices"]]$read(),
    p = g[["indptr"]]$read(),
    x = g[["data"]]$read(),
    dims = shape,
    index1 = FALSE
  )
  mat <- as(mat, "CsparseMatrix")

  gene_id <- g[["genes"]]$read()
  gene_symbol <- g[["gene_names"]]$read()
  barcodes <- g[["barcodes"]]$read()

  rownames(mat) <- gene_id
  colnames(mat) <- barcodes

  list(
    sample_id = sample_id,
    matrix = mat,
    features = data.frame(gene_id = gene_id, gene_symbol = gene_symbol, stringsAsFactors = FALSE),
    n_cells = ncol(mat)
  )
}

#' Loads every cohort C specimen listed in the sample manifest.
load_cohort_c_all <- function(sample_manifest, progress = TRUE) {
  c_manifest <- sample_manifest[sample_manifest$cohort == "C", , drop = FALSE]
  results <- vector("list", nrow(c_manifest))
  names(results) <- c_manifest$sample_id

  for (i in seq_len(nrow(c_manifest))) {
    sample_id <- c_manifest$sample_id[i]
    h5_path <- file.path(DATA_ROOT, c_manifest$file[i])
    if (progress) cat(sprintf("[C %d/%d] %s\n", i, nrow(c_manifest), sample_id))
    results[[sample_id]] <- load_cohort_c_sample(sample_id, h5_path)
  }
  results
}
