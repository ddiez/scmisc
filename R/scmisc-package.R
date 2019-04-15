#' scmisc.
#'
#' @name scmisc
#' @docType package
#'
#' @import tibble dplyr ggplot2 SummarizedExperiment SingleCellExperiment
#' @rawNamespace import(Biobase, except = combine)
#' @importFrom monocle newCellDataSet clusterCells reducedDimA<-
#' @importFrom Seurat CreateSeuratObject NormalizeData FindClusters GetAssayData Embeddings FetchData
#' @importFrom VGAM negbinomial.size
#' @importFrom MASS kde2d
#' @importFrom tidyr gather
#' @importFrom stats prcomp kmeans dist hclust cutree
#' @importFrom Rtsne Rtsne
#' @importFrom scran buildSNNGraph
#' @importFrom igraph cluster_louvain
NULL
