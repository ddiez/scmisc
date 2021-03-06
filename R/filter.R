#' filter_genes
#'
#' @param x object to filter.
#' @param min.cells minimum number of cells with min.count for gene to be preserved.
#' @param min.count minimum number of counts to be considered expressed.
#' @param assay.name name of assay to obtain matrix of counts.
#' @param ... arguments passed down to methods.
#'
#' @export
filter_genes <- function(x, ...) {
  UseMethod("filter_genes")
}

#' @rdname filter_genes
#' @export
filter_genes.SingleCellExperiment <- function(x, min.cells = 10, min.count = 0, assay.name = "counts", ...) {
  y <- assay(x, assay.name)
  sel.row <- rowSums(y > min.count) >= min.cells
  x[sel.row, ]
}

#' @rdname filter_genes
#' @export
filter_genes.matrix <- function(x, min.cells = 10, min.count = 0, ...) {
  sel.row <- rowSums(x > min.count) >= min.cells
  x[sel.row, ]
}

#' filter_cells
#'
#' @param x object to filter.
#' @param min.genes minimum number of genes with min.count for cell to be preserved.
#' @param min.count minimum number of counts to be considered expressed.
#' @param assay.name name of assay to obtain matrix of counts.
#' @param ... arguments passed down to methods.
#'
#' @export
filter_cells <- function(x, ...) {
  UseMethod("filter_cells")
}

#' @rdname filter_cells
#' @export
filter_cells.SingleCellExperiment <- function(x, min.genes = 10, min.count = 0, assay.name = "counts", ...) {
  y <- assay(x, assay.name)
  sel.col <- colSums(y > min.count) >= min.genes
  x[, sel.col]
}

#' @rdname filter_cells
#' @export
filter_cells.matrix <- function(x, min.genes = 10, min.count = 0, ...) {
  sel.col <- colSums(x > min.count) >= min.genes
  x[, sel.col]
}

