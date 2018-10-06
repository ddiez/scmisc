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
  d <- reducedDim(x, name) %>% fix_coords()

  if (annotate) {
    cdata <- colData(x)
    if (! is.null(cols))
      cdata <- cdata[, colnames(cdata) %in% cols, drop = FALSE]
    d <- cbind(d, cdata %>% as.data.frame())
  }
  d
}

#' @rdname get_coord
#' @export
get_coord.CellDataSet <- function(x, name = "A", annotate = TRUE, cols = NULL) {
  coord <- do.call(paste0("reducedDim", name), list(cds = x))
  d <- t(coord) %>% fix_coords()

  if (annotate) {
    pdata <- pData(x)
    if (! is.null(cols))
      pdata <- pdata[, colnames(pdata) %in% cols, drop = FALSE]
    d <- cbind(d, pdata %>% as.data.frame())
  }
  d
}

fix_coords <- function(x) {
  as.data.frame(x) %>% dplyr::rename(dim1 = 1, dim2 = 2)
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
  method <- match.arg(method, c("PCA", "tSNE"))

  y <- assay(x, assay.name)
  z <- reduce_dim(y, method = method, dims = dims, perplexity = perplexity, initial_dims = initial_dims, ...)
  reducedDim(x, coord.name) <- z

  x
}

#' @rdname reduce_dim
#' @export
reduce_dim.matrix <- function(x, method = "PCA", dims = 2, perplexity = NULL, initial_dims = 50, ...) {
  method <- match.arg(method, c("PCA", "tSNE"))

  if (method == "PCA") {
    dims <- seq_len(dims)
    z <- prcomp(t(x))[["x"]][, dims, drop = FALSE]
  }

  if (method == "tSNE") {
    if (is.null(perplexity)) {
      perplexity <- sqrt(ncol(x))
    }
    z <- Rtsne(t(x), dims = dims, perplexity = perplexity, initial_dims = initial_dims, ...)[["Y"]]
  }

  colnames(z) <- paste0("dim", seq_len(ncol(z)))
  rownames(z) <- colnames(x)

  z
}
