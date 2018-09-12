#' get_coord
#'
#' Extract coordinates (e.g. from dimensionality reduction) from bio datasets,
#' optionally with annotations to help plotting them with ggplot2 package.
#'
#' @param x object to obtain coordinates.
#' @param name name of the holder for coordinates.
#' @param annotate logical; whether to annotate each coordinate.
#' @param cols if not NULL the name of columns to include in annotation.
#'
#' @return A data.frame object.
#'
#' @export
get_coord <- function(x, name = "tSNE", annotate = TRUE, cols = NULL) {
  UseMethod("get_coord")
}

#' @rdname get_coord
#' @export
get_coord.SingleCellExperiment <- function(x, name = "tSNE", annotate = TRUE, cols = NULL) {
  d <- reducedDim(x, name) %>% as.data.frame() %>% dplyr::rename(dim1 = 1, dim2 = 2)

  if (annotate) {
    cdata <- colData(x)
    if (! is.null(cols))
      cdata <- cdata[, colnames(cdata) %in% cols, drop = FALSE]
    d <- cbind(d, cdata %>% as.data.frame())
  }
  d
}
