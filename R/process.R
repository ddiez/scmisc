#' process
#'
#' Applies several processing steps to a single cell genomics object.
#'
#' @param x an object of class Seurat.
#' @param assay assay to use for processing.
#' @param dims PCA dimensions to use for UMAP and clustering.
#' @param algorithm algorithm to use for clustering.
#' @param resolution resolution to use for clustering.
#' @param verbose whether to output diagnostic information.
#'
#' @export
process <- function(x, ...) {
  UseMethod("process")
}

#' @rdname process
#' @export
process.Seurat <- function(x, assay = NULL, dims = 1:10, algorithm = 1, resolution = 0.6, verbose = FALSE) {
  if (!is.null(assay)) {
    old.assay <- SeuratObject::DefaultAssay(x)
    DefaultAssay(x) <- assay
  }
  x <- FindVariableFeatures(x, verbose = verbose)
  x <- ScaleData(x, verbose = verbose)
  x <- RunPCA(x, verbose = verbose)
  x <- RunUMAP(x, dims = dims, verbose = verbose)
  x <- FindNeighbors(x, dims = dims, verbose = verbose)
  x <- FindClusters(x, algorithm = algorithm, resolution = resolution, verbose = verbose)

  if (!is.null(assay)) {
    DefaultAssay(x) <- old.assay
  }
  x
}
