#' cluster_cells
#'
#' @param x object with data matrix.
#' @param method type of dimensionality reduction.
#' @param ncluster number of clusters.
#' @param assay.name name of assay slot.
#' @param coord.name name of reducedDim slot.
#' @param column.name name of column to store cluster class.
#' @param ... arguments passed down to methods.
#'
#' @export
#'
cluster_cells <- function(x, method = "kmeans", ...) {
  UseMethod("cluster_cells")
}

#' @rdname cluster_cells
#' @export
cluster_cells.SingleCellExperiment <- function(x, method = "kmeans", ncluster = NULL, assay.name = "logcounts", coord.name = "PCA", column.name = paste0("cluster_", method), ...) {
  method <- match.arg(method, c("kmeans", "density"))

  y <- assay(x, assay.name)

  if (method == "kmeans") {
    cl <- kmeans(t(y), centers = ncluster, ...)
    colData(x)[[column.name]] <- factor(cl$cluster)
  }

  if (method == "density") {
    cds <- as.CellDataSet(x)
    reducedDimA(cds) <- t(reducedDim(x, coord.name))

    cds <- clusterCells(cds, num_clusters = ncluster, ...)
    colData(x)[[column.name]] <- factor(cds[["Cluster"]])
  }

  x
}
