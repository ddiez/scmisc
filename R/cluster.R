#' cluster_cells
#'
#' @param x object with data matrix.
#' @param method type of dimensionality reduction.
#' @param ncluster number of clusters.
#' @param assay.name name of assay slot.
#' @param coord.name name of reducedDim slot.
#' @param column.name name of column to store cluster class.
#' @param hclust.method method for hierarchical clustering (default: ward.D).
#' @param dist.method method for distance calculation (default: euclidean).
#' @param resolution parameter for seurat method (see ?Seurat::FindClusters).
#' @param algorithm parameter for seurat method (see ?Seurat::FindClusters).
#' @param ... arguments passed down to methods.
#'
#' @export
#'
cluster_cells <- function(x, method = "kmeans", ...) {
  UseMethod("cluster_cells")
}

#' @rdname cluster_cells
#' @export
cluster_cells.SingleCellExperiment <- function(x, method = "kmeans", ncluster = NULL, assay.name = "logcounts", coord.name = "PCA", column.name = "cluster", hclust.method = "ward.D", dist.method = "euclidean", resolution = 1, algorithm = 3, ...) {
  method <- match.arg(method, c("kmeans", "hclust", "louvain", "density", "leiden"))

  y <- assay(x, assay.name)

  if (method == "kmeans") {
    cl <- kmeans(t(y), centers = ncluster, ...)
    colData(x)[[column.name]] <- factor(cl[["cluster"]])
  }

  if (method == "hclust") {
    cl <- hclust(dist(t(y), method = dist.method), method = hclust.method, ...)
    colData(x)[[column.name]] <- factor(cutree(cl, k = ncluster))
  }

  if (method == "louvain") {
    g <- scran::buildSNNGraph(y)
    cl <- igraph::cluster_louvain(g)
    colData(x)[[column.name]] <- factor(cl$membership)
  }

  if (method == "density") {
    cds <- as.CellDataSet(x)
    monocle::reducedDimA(cds) <- t(reducedDim(x, coord.name))

    cds <- scran::clusterCells(cds, num_clusters = ncluster, ...)
    colData(x)[[column.name]] <- factor(cds[["Cluster"]])
  }

  if (method == "leiden") {
    g <- scran::buildSNNGraph(y)
    cl <- leiden::leiden(g, resolution_parameter = resolution)
    colData(x)[[column.name]] <- factor(cl)
  }

  x
}
