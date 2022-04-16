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

#' cluster_genes
#'
#' Cluster genes into modules using louvain/leiden algorithm
#'
#' @param x object with genes to clusters.
#' @param features which features to cluster.
#' @param reduction what reduction embeddings to use.
#' @param assay assay for Seurat objects.
#' @param slot slot for Seurat objects.
#' @param dims dimensions to use.
#' @param approx whether to use approximate method for PCA.
#' @param algorithm algorithm for clustering.
#' @param resolution resolution for clustering.
#' @param return.seurat whether to return a Seurat object.
#' @param ... arguments passed down to methods.
#'
#' @return
#' @export
#'
cluster_genes <- function(x, ...) {
  UseMethod("cluster_genes")
}

#' @rdname cluster_genes
#' @export
cluster_genes.Seurat <- function(x, features=NULL, reduction="pca", assay=NULL, slot="data", dims=1:10, approx=TRUE, algorithm="leiden", resolution=.6, return.seurat=FALSE, ...) {
  y <- SeuratObject::GetAssayData(x, assay=assay, slot=slot)
  if (!is.null(features))
    y <- y[features, , drop=FALSE]
  y <- Matrix::t(y)
  y <- SeuratObject::CreateSeuratObject(y)
  y <- Seurat::ScaleData(y, features=rownames(y), verbose=FALSE)
  y <- Seurat::RunPCA(y, features=rownames(y), approx=approx, verbose=FALSE)

  if (reduction=="umap") {
    y <- SeuratObject::RunUMAP(y, dims=dims, verbose=FALSE, metric="correlation")
    dims <- 1:2
  }

  y <- Seurat::FindNeighbors(y, reduction=reduction, dims=dims, verbose=FALSE)
  y <- Seurat::FindClusters(y, algorithm=algorithm, resolution=resolution, verbose=FALSE)

  if (return.seurat) {
    return(y)
  }
  else {
    y <- y[[]] |> rownames_to_column("gene") |> select(gene, module=seurat_clusters) |>
      mutate(module=fct_inseq(module))
  }
  y
}

