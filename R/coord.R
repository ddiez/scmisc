#' get_coord
#'
#' Extract coordinates (e.g. from dimensionality reduction) from bio datasets,
#' optionally with annotations to help plotting them with ggplot2 package.
#'
#' @param x object to obtain coordinates.
#' @param coord.name name of the holder for coordinates.
#' @param add.cols logical; whether to annotate each coordinate with column data (i.e. colData).
#' @param add.exprs if not NULL a list of genes to add counts obtained with get_expression().
#'
#' @return A data.frame object.
#'
#' @export
get_coord <- function(x, coord.name = "TSNE", add.cols = TRUE, add.exprs = NULL) {
  UseMethod("get_coord")
}

#' @rdname get_coord
#' @export
get_coord.SingleCellExperiment <- function(x, coord.name = "TSNE", add.cols = TRUE, add.exprs = FALSE) {
  d <- reducedDim(x, coord.name)[, 1:2] %>% fix_coords()

  if (! isFALSE(add.cols)) {
    cdata <- colData(x)
    if (! isTRUE(add.cols))
      cdata <- cdata[, colnames(cdata) %in% add.cols, drop = FALSE]
    d <- cbind(d, cdata %>% as.data.frame())
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
get_coord.CellDataSet <- function(x, coord.name = "A", add.cols = TRUE, add.exprs = NULL) {
  coord <- do.call(paste0("reducedDim", coord.name), list(cds = x))
  d <- t(coord)[, 1:2] %>% fix_coords()

  if (! isFALSE(add.cols)) {
    pdata <- pData(x)
    if (! isTRUE(add.cols))
      pdata <- pdata[, colnames(pdata) %in% add.cols, drop = FALSE]
    d <- cbind(d, pdata %>% as.data.frame())
  }
  d
}

#' @rdname get_coord
#' @export
get_coord.seurat <- function(x, coord.name = "tsne", add.cols = TRUE, add.exprs = NULL) {
  d <- x@dr[[coord.name]]@cell.embeddings[, 1:2] %>% fix_coords()

  if (! isFALSE(add.cols)) {
    cdata <- x@meta.data
    if (! isTRUE(add.cols))
      cdata <- cdata[, colnames(cdata) %in% add.cols, drop = FALSE]
    d <- cbind(d, cdata)
  }
  d
}

fix_coords <- function(x) {
  as.data.frame(x) %>% dplyr::rename(dim1 = 1, dim2 = 2)
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
  if (! col.name %in% colnames(x)) stop("col.name not a colname of x.")

  ll <- unique(x[[col.name]])
  tmp <- sapply(ll, function(l) {
    x[[col.name]] == l
  }) %>% as.data.frame()
  colnames(tmp) <- ll
  y <- bind_cols(tmp, x)

  y <- y %>% gather(!!out.name, !!value.name, seq_len(length(ll)))

  if (is.factor(x[[col.name]]))
    y[[out.name]] <- factor(y[[out.name]], levels = levels(x[[col.name]]))

  if (is.numeric(x[[col.name]]))
    y[[out.name]] <- as.numeric(y[[out.name]])

  if (is.integer(x[[col.name]]))
    y[[out.name]] <- as.integer(y[[out.name]])

  y
}

#' reduce_dim
#'
#' @param x matrix object.
#' @param method tyoe of dimensionality reduction.
#' @param dims number of dimensions to keep in final result.
#' @param assay.name name of assay slot.
#' @param coord.name name of reducedDim slot.
#' @param perplexity perplexity for tSNE.
#' @param initial_dims initial PCA dimensions for tSNE.
#' @param ...  arguments to be passed down to methods.
#'
#' @export
#'
reduce_dim <- function(x, method = "PCA", dims = 2, ...) {
  UseMethod("reduce_dim")
}

#' @rdname reduce_dim
#' @export
reduce_dim.SingleCellExperiment <- function(x, method = "PCA", dims = 2, assay.name = "logcounts", coord.name = method, perplexity = NULL, initial_dims = 50, ...) {
  method <- match.arg(method, c("PCA", "TSNE"))

  y <- assay(x, assay.name)
  z <- reduce_dim(y, method = method, dims = dims, perplexity = perplexity, initial_dims = initial_dims, ...)
  reducedDim(x, coord.name) <- z

  x
}

#' @rdname reduce_dim
#' @export
reduce_dim.matrix <- function(x, method = "PCA", dims = 2, perplexity = NULL, initial_dims = 50, ...) {
  method <- match.arg(method, c("PCA", "TSNE"))

  if (method == "PCA") {
    dims <- seq_len(dims)
    z <- prcomp(t(x))[["x"]][, dims, drop = FALSE]
  }

  if (method == "TSNE") {
    if (is.null(perplexity)) {
      perplexity <- sqrt(ncol(x))
    }
    z <- Rtsne(t(x), dims = dims, perplexity = perplexity, initial_dims = initial_dims, ...)[["Y"]]
  }

  colnames(z) <- paste0("dim", seq_len(ncol(z)))
  rownames(z) <- colnames(x)

  z
}
