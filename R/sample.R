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
sample_cells.Seurat <- function(x, group = NULL, n = 20, ...) {
  cdata <- x@meta.data %>%
    as_tibble(rownames = ".id")

  if (!is.null(group))
    cdata <- cdata %>% group_by_(group)

  sel.cells <- cdata %>%
    sample_n(n) %>%
    pull(.data$.id)

  x[, sel.cells]
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


#' create_training_sets
#'
#' Divides a SingleCellExperiment object into train/test subsets.
#'
#' @param x SingleCellExperiment object.
#' @param frac fraction of dataset for testing.
#'
#' @return list with two SingleCellExperiment, one with train and other with test data.
#' @export
create_training_sets <- function(x, frac = .1) {
  cdata <- colData(x) %>% as.data.frame()
  id.test <- cdata %>% sample_frac(frac) %>% pull("cell_index")
  id.train <- cdata %>% filter(! .data$cell_index %in% id.test) %>% pull("cell_index")
  list(train = x[, id.train], test = x[, id.test])
}
