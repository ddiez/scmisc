#' sort_features
#'
#' Sort features (typically genes) using Pearson correlation.
#'
#' @param x an object with features to sort.
#' @param features a character vector of features to sort.
#' @param method method to be used for sortening.
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
sort_features.Seurat <- function(x, features = NULL, method = "pearson", slot = "data", assay = NULL, ...) {
  if (is.null(assay)) assay <- DefaultAssay(x)

  m <- GetAssayData(x, slot = slot, assay = assay)

  if (!is.null(features)) {
    features <- features[features %in% rownames(m)]
    m <- m[features, ]
  }
  sort_features(as.matrix(m), method = method)
}

#' @rdname sort_features
#' @export
sort_features.SingleCellExperiment <- function(x, features = NULL, method = "pearson", assay = NULL, ...) {
  if (is.null(assay)) assay <- names(assays(x))[1]

  m <- assay(x, assay)

  if (!is.null(features)) {
    features <- features[features %in% rownames(m)]
    m <- m[features, ]
  }
  sort_features(m, method = method)
}

#' @rdname sort_features
#' @export
sort_features.matrix <- function(x, method = "pearson", ...) {
  method <- match.arg(method, c("pearson", "sum"))
  switch(method,
         "pearson" = {
           ord <- hclust(as.dist(t(1 - cor(t(x)))))$order
           features <- rownames(x)[ord]
         },
         "sum" = {
           features <- names(sort(rowSums(x)))
         }
  )
  features
}
