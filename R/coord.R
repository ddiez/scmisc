#' get_coord
#'
#' Extract coordinates (e.g. from dimensionality reduction) from bio datasets,
#' optionally with annotations to help plotting them with ggplot2 package.
#'
#' @param x object to obtain coordinates.
#' @param coord.name name of the holder for coordinates.
#' @param add.cols logical; whether to annotate each coordinate with column data (i.e. colData).
#' @param add.exprs if not NULL a list of genes to add counts obtained with get_expression().
#' @param assay name of assay for Seurat objects.
#' @param slot name of slot for Seurat objects.
#' @param ... arguments passed down to methods.
#'
#' @return A data.frame object.
#'
#' @export
get_coord <- function(x, coord.name = NULL, add.cols = TRUE, add.exprs = NULL, ...) {
  UseMethod("get_coord")
}

#' @rdname get_coord
#' @export
get_coord.SingleCellExperiment <- function(x, coord.name = NULL, add.cols = TRUE, add.exprs = FALSE, ...) {
  if (is.null(coord.name)) {
    coord.name <- names(reducedDims(x))[1]
  }

  if (is.null(coord.name)) {
    stop("No coordinates found in object. Run reduce_dim() to generate some.")
  }

  d <- reducedDim(x, coord.name)[, 1:2] |> fix_coords()

  if (! isFALSE(add.cols)) {
    cdata <- colData(x)
    if (! isTRUE(add.cols))
      cdata <- cdata[, colnames(cdata) %in% add.cols, drop = FALSE]
    d <- cbind(d, cdata |> as.data.frame())
  }

  if (! isFALSE(add.exprs)) {
    if (isTRUE(add.exprs)) {
      add.exprs <- rowData(x)["symbol"]
    }
    d <- cbind(d, sapply(add.exprs, get_expression, x = x))
  }

  d
}

#' @rdname get_coord
#' @export
get_coord.CellDataSet <- function(x, coord.name = "A", add.cols = TRUE, add.exprs = FALSE, ...) {
  coord <- do.call(paste0("reducedDim", coord.name), list(cds = x))
  d <- t(coord)[, 1:2] |> fix_coords()

  if (! isFALSE(add.cols)) {
    pdata <- pData(x)
    if (! isTRUE(add.cols))
      pdata <- pdata[, colnames(pdata) %in% add.cols, drop = FALSE]
    d <- cbind(d, pdata |> as.data.frame())
  }
  d
}

#' @rdname get_coord
#' @export
get_coord.seurat <- function(x, coord.name = "tsne", add.cols = TRUE, add.exprs = FALSE, ...) {
  d <- x@dr[[coord.name]]@cell.embeddings[, 1:2] |> fix_coords()

  if (! isFALSE(add.cols)) {
    cdata <- x@meta.data
    if (! isTRUE(add.cols))
      cdata <- cdata[, colnames(cdata) %in% add.cols, drop = FALSE]
    d <- cbind(d, cdata)
  }
  d
}

#' @rdname get_coord
#' @export
get_coord.Seurat <- function(x, coord.name = NULL, add.cols = TRUE, add.exprs = FALSE, assay = NULL, slot = "data", ...) {
  if (is.null(coord.name)) {
    coord.name <- Reductions(x)[1]
  }

  if (length(coord.name) == 0) {
    stop("No coordinates found in object. Run reduce_dim() to generate some.")
  }

  d <- Embeddings(x, coord.name)[, 1:2] |> fix_coords()

  if (! isFALSE(add.cols)) {
    cdata <- x@meta.data
    if (! isTRUE(add.cols))
      cdata <- cdata[, colnames(cdata) %in% add.cols, drop = FALSE]
    d <- cbind(d, cdata)
  }

  if (! isFALSE(add.exprs)) {
    if (isTRUE(add.exprs)) {
      add.exprs <- rownames(GetAssayData(x, assay = assay, slot = slot))
    }
    d <- cbind(d, sapply(add.exprs, get_expression, x = x, assay = assay, slot = slot))
  }

  d
}

fix_coords <- function(x) {
  as.data.frame(x) |> dplyr::rename(dim1 = 1, dim2 = 2)
}

#' expand_column
#'
#' @param x a data.frame or similar object.
#' @param col.name  column from the data.frame to expand.
#' @param out.name  column in the data.frame where to store expanded values.
#' @param value.name  column in the data.frame where to store logical vector.
#' @param ... arguments passed down to methods.
#'
#' @export
expand_column <- function(x, ...) {
  UseMethod("expand_column")
}

#' @rdname expand_column
#' @export
expand_column.data.frame <- function(x, col.name = NULL, out.name = col.name, value.name = "value", ...) {
  if (is.null(col.name)) stop("col.name is required.")
  if (! all(col.name %in% colnames(x))) stop("col.name not a colname of x.")

  if (length(col.name) > 2) stop("1 or 2 column names are required.")

  if (length(col.name) == 1)
    d <- expand_column1d(x, col.name, out.name = col.name, value.name = value.name)
  if (length(col.name) == 2)
    d <- expand_column2d(x, col.name)

  d
}


expand_column1d <- function(x, col.name = NULL, out.name = col.name, value.name = "value") {
  if (length(col.name) != 1) stop("col.names must be a character vector of length 2.")

  ll <- unique(x[[col.name]])
  tmp <- sapply(ll, function(l) {
    x[[col.name]] == l
  }) |> as.data.frame()
  colnames(tmp) <- ll
  y <- bind_cols(tmp, x)

  y <- y |> gather(!!out.name, !!value.name, seq_len(length(ll)))

  if (is.factor(x[[col.name]]))
    y[[out.name]] <- factor(y[[out.name]], levels = levels(x[[col.name]]))

  if (is.numeric(x[[col.name]]))
    y[[out.name]] <- as.numeric(y[[out.name]])

  if (is.integer(x[[col.name]]))
    y[[out.name]] <- as.integer(y[[out.name]])

  y
}

expand_column2d <- function(x, col.names) {
  if (length(col.names) != 2) stop("col.names must be a character vector of length 2.")
  val.cols <- paste0("val.", col.names)
  names(val.cols) <- col.names
  d <- x
  for (col in col.names) {
    d <- expand_column1d(d, col, value.name = val.cols[col])
  }
  d |> mutate(value = .data[[val.cols[1]]] & .data[[val.cols[2]]]) |>
    arrange(.data[["value"]])
}

#' reduce_dim
#'
#' @param x matrix object.
#' @param method tyoe of dimensionality reduction: pca, tsne, umap.
#' @param dims number of dimensions to keep in final result.
#' @param assay.name name of assay slot.
#' @param coord.name name of reducedDim slot.
#' @param perplexity perplexity for tSNE.
#' @param initial_dims initial PCA dimensions for tSNE.
#' @param ...  arguments to be passed down to methods.
#'
#' @export
#'
reduce_dim <- function(x, method = "pca", dims = 2, ...) {
  UseMethod("reduce_dim")
}

#' @rdname reduce_dim
#' @export
reduce_dim.SingleCellExperiment <- function(x, method = "pca", dims = 2, assay.name = "logcounts", coord.name = method, perplexity = NULL, initial_dims = 50, ...) {
  method <- match.arg(method, c("pca", "tsne", "umap"))

  y <- assay(x, assay.name)
  z <- reduce_dim(y, method = method, dims = dims, perplexity = perplexity, initial_dims = initial_dims, ...)
  reducedDim(x, coord.name) <- z

  x
}

#' @rdname reduce_dim
#' @export
reduce_dim.matrix <- function(x, method = "pca", dims = 2, perplexity = 30, initial_dims = 50, ...) {
  method <- match.arg(method, c("pca", "tsne", "umap"))

  if (method == "pca") {
    dims <- seq_len(dims)
    z <- prcomp(t(x))[["x"]][, dims, drop = FALSE]
  }

  if (method == "tsne") {
    z <- Rtsne::Rtsne(t(x), dims = dims, perplexity = perplexity, initial_dims = initial_dims, ...)[["Y"]]
  }

  if (method == "umap") {
    z <- uwot::umap(t(x), ret_model = FALSE, ...)
  }

  colnames(z) <- paste0("dim", seq_len(ncol(z)))
  rownames(z) <- colnames(x)

  z
}

#' get_coord_names
#'
#' @param x object with coordinates
#'
#' @export
get_coord_names <- function(x) {
  UseMethod("get_coord_names")
}


#' @rdname get_coord_names
#' @export
get_coord_names.SingleCellExperiment <- function(x) {
  reducedDimNames(x)
}

#' @rdname get_coord_names
#' @export
get_coord_names.seurat <- function(x) {
  names(x$dr)
}
