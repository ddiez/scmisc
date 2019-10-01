#' scmisc.
#'
#' @name scmisc
#' @docType package
#'
#' @import tibble dplyr ggplot2 ggraph tidygraph SummarizedExperiment SingleCellExperiment
#' @rawNamespace import(Biobase, except = combine)
#' @importFrom monocle newCellDataSet clusterCells reducedDimA<-
#' @importFrom Seurat DefaultAssay GetAssayData Embeddings FetchData Reductions
#' @importFrom slingshot SlingshotDataSet slingPseudotime slingParams slingAdjacency
#' @importFrom VGAM negbinomial.size
#' @importFrom MASS kde2d
#' @importFrom tidyr gather
#' @importFrom stats prcomp kmeans dist hclust cutree
#' @importFrom Rtsne Rtsne
#' @importFrom uwot umap
#' @importFrom leiden leiden
#' @importFrom scran buildSNNGraph
#' @importFrom igraph cluster_louvain graph_from_adjacency_matrix
#' @importFrom ggbeeswarm geom_quasirandom
NULL
