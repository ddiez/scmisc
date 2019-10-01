#' sort_features
#'
#' Sort features (typically genes) using Pearson correlation.
#'
#' @param x an object with features to sort.
#' @param features a character vector of features to sort.
#' @param slot the slot to obtain the data from.
#' @param assay the assay to obtain the data from.
#' @param ... rguments passed down to methods.
#'
#' @export
sort_features <- function(x, ...) {
  UseMethod("sort_features")
}

#' @rdname sort_features
#' @export
sort_features.Seurat <- function(x, features = NULL, slot = "scale.data", assay = NULL, ...) {
  if (is.null(assay)) assay <- DefaultAssay(x)

  m <- GetAssayData(x, slot = slot, assay = assay)

  if (!is.null(features)) {
    features <- features[features %in% rownames(m)]
    m <- m[features, ]
  }
  sort_features(m)
}

#' @rdname sort_features
#' @export
sort_features.SingleCellExperiment <- function(x, features = NULL, assay = NULL, ...) {
  if (is.null(assay)) assay <- names(assays(x))[1]

  m <- assay(x, assay)

  if (!is.null(features)) {
    features <- features[features %in% rownames(m)]
    m <- m[features, ]
  }
  sort_features(m)
}

#' @rdname sort_features
#' @export
sort_features.matrix <- function(x, ...) {
  ord <- hclust(as.dist(t(1 - cor(t(x)))))$order
  rownames(x)[ord]
}
