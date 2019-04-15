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
    as_tibble(rownames = ".id")

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
#' @param ... arguments passed down to methods (currently unused).
#'
#' @return list with two SingleCellExperiment, one with train and other with test data.
#' @export
create_training_sets <- function(x, ...) {
  UseMethod("create_training_sets")
}

#' @rdname create_training_sets
#' @export
create_training_sets.SingleCellExperiment <- function(x, frac = .1, ...) {
  cdata <- colData(x) %>% as_tibble(rownames = ".id")
  id.test <- cdata %>% sample_frac(frac) %>% pull(".id")
  id.train <- cdata %>% filter(! .data[[".id"]] %in% id.test) %>% pull("cell_index")
  list(train = x[, id.train], test = x[, id.test])
}

#' @rdname create_training_sets
#' @export
create_training_sets.Seurat <- function(x, frac = .1, ...) {
  cdata <- x@meta.data %>% as_tibble(rownames = ".id")
  id.test <- cdata %>% sample_frac(frac) %>% pull(".id")
  id.train <- cdata %>% filter(! .data[[".id"]] %in% id.test) %>% pull(".id")
  list(train = x[, id.train], test = x[, id.test])
}
