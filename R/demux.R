#' demux
#'
#' @param x object to demux.
#' @param n number of hashtags.
#' @param ... arguments passed down to methods.
#'
#' @export
demux <- function(x, n, ...) {
  UseMethod("demux")
}

#' @rdname demux
#' @export
demux.Seurat <- function(x, n, assay = NULL) {
  if (is.null(assay))
    assay <- DefaultAssay(x)

  m <- GetAssayData(x, assay = assay, slot = "counts")
  m <- log2(m + 1)
  m <- t(scale(m))
  m <- as.data.frame(m)
  m <- na.omit(m)


  ms <- t(apply(m, 1, function(xx) sort(xx, decreasing = TRUE)[1:2]))
  colnames(ms) <- c("first", "second")

  m <- cbind(m, ms)
  m$diff <- m$first - m$second
  m <- m[m$diff > 0.1, ]
  res <- cluster::pam(m, k = n)
  m$cluster <- as.character(res$clustering)
  m
}

#' demux_cutoff
#'
#' Hashtag demultiplexing using cutoffs, decided with plot_pairs().
#'
#' @param x Seurat object.
#' @param cutoff vector of cutoffs to apply to each hashtag.
#' @param assay assay from Seurat object.
#' @param slot slot from Seurat object.
#'
#' @export
#'
demux_cutoff <- function(x, cutoff, assay = "HTO", slot = "data") {
  d <- t(GetAssayData(x, assay = assay, slot = slot))
  d <- t(t(d) > cutoff)

  apply(d, 1, function(x) {
    if (sum(x) == 0)
      return("Negative")

    if (sum(x) > 1)
      return("Doublet")

    return(colnames(d)[x])
  })
}
