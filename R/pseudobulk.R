#' pseudobulk
#'
#' Generate pseudobulk counts of single cell data per cell population and sample name. Return a list of DGEList objects.
#'
#' @param x Seurat object.
#' @param split.by Split single cell data using this column. Typically cell populations or clusters.
#' @param group.by Aggregate counts from cells in each group. Typically sample names.
#' @param samples A data.frame with sample information to include in slot "samples" in the DGEList object.
#' @param genes A data.frame with gene information to include in slot "genes" in the DGEList object.
#' @param assay Seurat assay to use. Default: NULL (default assay).
#' @param layers Seurat layers to use. Default: counts.
#' @param sort_samples whetner to sort the columns of the pseudobulk samples based on the order in the samples data.frame.
#'
#' @export
pseudobulk <- function(x, split.by, group.by, samples=NULL, genes=NULL, sort_samples=TRUE, ...) {
  UseMethod("pseudobulk")
}

#' @rdname pseudobulk
#' @export
pseudobulk.Seurat <- function(x, split.by, group.by, samples=NULL, genes=NULL, assay=NULL, layers="counts", sort_samples=FALSE, ...) {
  if (!is.null(assay))
    DefaultAssay(x) <- assay

  assay <- DefaultAssay(x)

  if (packageVersion("Seurat") >= "5.0.0")
    x <- Seurat::DietSeurat(x, layers=layers, assay=assay)
  else
    x <- Seurat::DietSeurat(x, counts=TRUE, data=TRUE, assay=assay)

  groups <- x[[]][[split.by]]
  if (!is.factor(groups)) groups <- factor(groups)
  group_levels <- levels(groups)
  xl <- lapply(group_levels, function(group) {
    x[, groups == group]
  })
  names(xl) <- group_levels

  if (is.null(samples)) {
    samples <- x[[]][[group.by]]
    if (!is.factor(samples))
      samples <- factor(samples)

    samples <- levels(samples)
    samples <- data.frame(row.names=samples)
  }

  xl <- lapply(xl, function(x) {
    m <- matrix(0, nrow=nrow(x), ncol=nrow(samples), dimnames=list(rownames(x), rownames(samples)))

    if (packageVersion("Seurat") >= "5.0.0")
      tmp <- Seurat::AggregateExpression(x, group.by=group.by, assay=assay)[[assay]]
    else
      tmp <- Seurat:::PseudobulkExpression(x, pb.method="aggregate", group.by=group.by, assay=assay, slot=layers[1])[[assay]]

    m[rownames(tmp), colnames(tmp)] <- tmp
    m
  })

  dge <- lapply(names(xl), function(n) {
    tmp <- edgeR::DGEList(xl[[n]], samples=samples, genes=genes)
    tmp$samples$cluster <- n
    tmp
  })
  names(dge) <- names(xl)

  if (sort_samples) {
    dge <- lapply(dge, function(x) {
      x[, sort(rownames(samples))]
    })
  }

  dge
}
