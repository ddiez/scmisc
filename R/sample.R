#' sample_cells
#'
#' Subset a single cell experiment object by sampling the number of specified cells,
#' optionally by grouping.
#'
#' @param x an object to sample from.
#' @param group optional grouping variable.
#' @param n number of cells to sample (per group).
#' @param ... parameters passed down to methods.
#'
#' @export
sample_cells <- function(x, ...) {
  UseMethod("sample_cells")
}

#' @rdname sample_cells
#' @export
sample_cells.SingleCellExperiment <- function(x, group = NULL, n = 20, ...) {
  cdata <- colData(x) %>%
    as.data.frame() %>%
    as.tibble(rownames = ".id")

  if (!is.null(group))
    cdata <- cdata %>% group_by_(group)

  sel.cells <- cdata %>%
    sample_n(n) %>%
    pull(.data$.id)

  x[, sel.cells]
}
