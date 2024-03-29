#' sample_cells
#'
#' Subset a single cell experiment object by sampling the number of specified cells,
#' optionally by grouping. If neither n, frac or n_max is specified the lowest number
#' of cells in group will be used.
#'
#' @param x an object to sample from.
#' @param group optional grouping variable.
#' @param n number of cells to sample (per group).
#' @param frac fraction of cells to sample (per group).
#' @param n_max number of cells to sample, if less cells available then all cells returned.
#' @param ... parameters passed down to methods.
#'
#' @export
sample_cells <- function(x, ...) {
  UseMethod("sample_cells")
}

#' @rdname sample_cells
#' @export
sample_cells.Seurat <- function(x, group = NULL, n = NULL, frac = NULL, n_max = NULL, ...) {
  cdata <- x[[]] |>
    as_tibble(rownames = ".id")

  if (!is.null(group))
    cdata <- cdata |> group_by(.data[[group]])

  if (is.null(n) && is.null(frac) && is.null(n_max)) {
    n <- cdata |> count() |> pull(n) |> min()
    message("Neither 'n', 'frac' or 'n_max' specified. Sampling ", n, " cells per group.")
  }

  if (!is.null(n)) {
    sel.cells <- cdata |>
      sample_n(n) |>
      pull(.data[[".id"]])
  }

  if (!is.null(frac)) {
    sel.cells <- cdata |>
      sample_frac(frac) |>
      pull(.data[[".id"]])
  }

  if (!is.null(n_max)) {
    d <- split(cdata, cdata[[group]])
    sel.cells <- lapply(d, function(x) if (nrow(x) < n_max) x else sample_n(x, size = n_max)) |>
      bind_rows() |>
      pull(.data[[".id"]])
  }

  x[, sel.cells]
}

#' @rdname sample_cells
#' @export
sample_cells.SingleCellExperiment <- function(x, group = NULL, n = NULL, frac = NULL, ...) {
  cdata <- colData(x) |>
    as.data.frame() |>
    as_tibble(rownames = ".id")

  if (!is.null(group))
    cdata <- cdata |> group_by_at(.vars = group)

  if (is.null(n) && is.null(frac)) {
    n <- cdata |> count() |> pull(n) |> min()
    message("Neither 'n' or 'frac' specified. Sampling ", n, " cells per group.")
  }

  if (!is.null(n)) {
    sel.cells <- cdata |>
      sample_n(n) |>
      pull(.data$.id)
  }

  if (!is.null(frac)) {
    sel.cells <- cdata |>
      sample_frac(frac) |>
      pull(.data$.id)
  }

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
  cdata <- colData(x) |> as_tibble(rownames = ".id")
  id.test <- cdata |> sample_frac(frac) |> pull(".id")
  id.train <- cdata |> filter(! .data[[".id"]] %in% id.test) |> pull("cell_index")
  list(train = x[, id.train], test = x[, id.test])
}

#' @rdname create_training_sets
#' @export
create_training_sets.Seurat <- function(x, frac = .1, ...) {
  cdata <- x@meta.data |> as_tibble(rownames = ".id")
  id.test <- cdata |> sample_frac(frac) |> pull(".id")
  id.train <- cdata |> filter(! .data[[".id"]] %in% id.test) |> pull(".id")
  list(train = x[, id.train], test = x[, id.test])
}
