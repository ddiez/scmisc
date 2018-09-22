#' Compute summary from assay matrix
#'
#' @param sum.fun summary function (by default mean).
#' @param name name of the feature value for which to compute the summary.
#' @param column name of the column where to find 'name' (default: 'symbol').
#' @param assay.name name of the assay slot (defualt: 'logcounts').
#' @param ... arguments passed down to specific methods.
#'
#' @export
get_assay_summary <- function(x, sum.fun = mean, ...) {
  UseMethod("get_assay_summary")
}

#' @rdname get_assay_summary
#' @export
get_assay_summary.SingleCellExperiment <- function(x, name, assay.name = "logcounts", column = "symbol", sum.fun = mean) {
  y <- assay(x, assay.name)
  z <- y[rowData(x)[[column]] == name, , drop = FALSE]
  apply(z, 2, sum.fun)
}

#' @rdname get_assay_summary
#' @export
get_expression <- function(...) {
  get_assay_summary(...)
}
