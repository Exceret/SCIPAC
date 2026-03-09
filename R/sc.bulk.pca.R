#' Use PCA to perform dimension reduction on single-cell and bulk data
#'
#' @param sc.dat a normalized, log-transformed, and hvg's-selected single-cell data matrix, whose rows are for genes and columns are for cells.
#' @param bulk.dat a normalized and log-transformed bulk data matrix, whose rows are for genes and columns are for cells.
#' This data matrix uses the same gene set as \code{sc.dat}.
#' @param do.pca.sc if \code{TRUE}, first do PCA on \code{sc.dat} and use the rotation matrix on \code{bulk.dat};
#' if \code{False}, first do PCA on \code{bulk.data} and use the rotation matrix on \code{sc.dat}. The default is \code{FALSE}.
#' @param n.pc the number of PCs used for downstream analysis. The default is \code{60}.
#' @param batch_var a vector indicating batch assignment for each cell in sc.dat. If provided, harmony batch correction will be applied. The default is \code{NULL}.
#'
#' @return A list object contains dimension-reduced single-cell and bulk data.
#' @export

sc.bulk.pca <- function(
  sc.dat,
  bulk.dat,
  do.pca.sc = FALSE,
  n.pc = 60,
  batch_var = NULL
) {
  # Check batch_var if provided
  if (!is.null(batch_var)) {
    if (length(batch_var) != ncol(sc.dat)) {
      stop("Length of batch_var must match the number of cells in sc.dat")
    }
    if (!requireNamespace("harmony", quietly = TRUE)) {
      stop(
        "Package 'harmony' is required for batch correction. Please install it using: install.packages('harmony')"
      )
    }
  }

  if (all(rownames(sc.dat) == rownames(bulk.dat))) {
    # Do the centering
    sc.dat.cen <- scale(
      t(sc.dat),
      center = rowSums(sc.dat) / ncol(sc.dat),
      scale = FALSE
    )
    bulk.dat.cen <- scale(
      t(bulk.dat),
      center = rowSums(bulk.dat) / ncol(bulk.dat),
      scale = FALSE
    )

    if (do.pca.sc) {
      # Do PCA from single cell data
      sc.dat.pca <- stats::prcomp(
        sc.dat.cen,
        center = FALSE,
        scale. = FALSE
      )
      sc.dat.pca.rotation <- sc.dat.pca[["rotation"]]

      sc.dat.rot <- sc.dat.cen %*% sc.dat.pca.rotation
      bulk.dat.rot <- bulk.dat.cen %*% sc.dat.pca.rotation

      # Apply harmony if batch_var is provided
      if (!is.null(batch_var)) {
        meta_data <- data.frame(
          batch = batch_var,
          row.names = rownames(sc.dat.rot)
        )
        sc.dat.rot.corrected <- harmony::RunHarmony(
          data_mat = sc.dat.rot[, 1:n.pc],
          meta_data = meta_data,
          vars_use = "batch"
        )
        return(list(
          "sc.dat.rot" = sc.dat.rot.corrected,
          "bulk.dat.rot" = bulk.dat.rot[, 1:n.pc]
        ))
      }

      return(list(
        "sc.dat.rot" = sc.dat.rot[, 1:n.pc],
        "bulk.dat.rot" = bulk.dat.rot[, 1:n.pc]
      ))
    } else {
      # Do PCA from bulk data
      bulk.dat.pca <- stats::prcomp(
        bulk.dat.cen,
        center = FALSE,
        scale. = FALSE
      )
      bulk.dat.pca.rotation <- bulk.dat.pca[["rotation"]]

      sc.dat.rot <- sc.dat.cen %*% bulk.dat.pca.rotation
      bulk.dat.rot <- bulk.dat.cen %*% bulk.dat.pca.rotation

      # Apply harmony if batch_var is provided
      if (!is.null(batch_var)) {
        meta_data <- data.frame(
          batch = batch_var,
          row.names = rownames(sc.dat.rot)
        )
        sc.dat.rot.corrected <- harmony::RunHarmony(
          data_mat = sc.dat.rot[, 1:n.pc],
          meta_data = meta_data,
          vars_use = "batch"
        )
        return(list(
          "sc.dat.rot" = sc.dat.rot.corrected,
          "bulk.dat.rot" = bulk.dat.rot[, 1:n.pc]
        ))
      }

      return(list(
        "sc.dat.rot" = sc.dat.rot[, 1:n.pc],
        "bulk.dat.rot" = bulk.dat.rot[, 1:n.pc]
      ))
    }
  } else {
    stop("Row names of sc.dat and bulk.data do not match.")
  }
}
