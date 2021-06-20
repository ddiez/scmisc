#' as.CellDataSet
#'
#' Convert several data types into monocle::CellDataSet object.
#'
#' @param x object to convert.
#' @param expressionFamily expression family function to use (default: negbinomial.size)
#' @param gene_short_name column to use to set gene_short_name
#'
#' @return a CellDataSet object.
#' @export
as.CellDataSet <- function(x, expressionFamily = NULL, gene_short_name = "symbol") {
  UseMethod("as.CellDataSet")
}

#' @rdname as.CellDataSet
#' @export
as.CellDataSet.SingleCellExperiment <- function(x, expressionFamily = NULL, gene_short_name = "symbol") {
  exprs <- assay(x, "counts")
  pdata <- AnnotatedDataFrame(as.data.frame(colData(x)))
  fdata <- AnnotatedDataFrame(as.data.frame(rowData(x)))
  rownames(fdata) <- rownames(exprs)

  fdata[["gene_short_name"]] <- fdata[[gene_short_name]]

  if (is.null(expressionFamily))
    expressionFamily <- VGAM::negbinomial.size()

  monocle::newCellDataSet(exprs, phenoData = pdata, featureData = fdata, expressionFamily = expressionFamily)
}
