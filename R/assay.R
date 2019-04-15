#' Compute summary from assay matrix
#'
#' @param x object to compute summary.
#' @param sum.fun summary function (by default mean).
#' @param name name of the feature value for which to compute the summary.
#' @param column name of the column where to find 'name' (default: 'symbol').
#' @param assay.name name of the assay slot.
#' @param slot name of the Seurat data slot.
#' @param ... arguments passed down to specific methods.
#'
#' @export
get_assay_summary <- function(x, ...) {
  UseMethod("get_assay_summary")
}

#' @rdname get_assay_summary
#' @export
get_assay_summary.SingleCellExperiment <- function(x, name, assay.name = "logcounts", column = "symbol", sum.fun = mean, ...) {
  y <- assay(x, assay.name)
  z <- y[rowData(x)[[column]] == name, , drop = FALSE]
  apply(z, 2, sum.fun)
}

#' @rdname get_assay_summary
#' @export
get_assay_summary.Seurat <- function(x, name, assay.name = "RNA", slot = "data", column = "symbol", sum.fun = mean, ...) {
  y <- GetAssayData(x, slot, assay.name)
  z <- y[get_rowdata(x)[[column]] == name, , drop = FALSE]
  apply(z, 2, sum.fun)
}

#' @rdname get_assay_summary
#' @export
get_expression <- function(...) {
  get_assay_summary(...)
}
