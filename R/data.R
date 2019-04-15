#' get_coldata
#'
#' @param x object.
#'
#' @export
get_coldata <- function(x) {
  UseMethod("get_coldata")
}

#' @rdname get_coldata
#' @export
get_coldata.SingleCellExperiment <- function(x) {
  colData(x) %>% as.data.frame()
}

#' @rdname get_coldata
#' @export
get_coldata.seurat <- function(x) {
  x@meta.data %>% as.data.frame()
}

#' @rdname get_coldata
#' @export
get_coldata.Seurat <- function(x) {
  x@meta.data %>% as.data.frame()
}


#' get_rowdata
#'
#' @param x object.
#'
#' @export
get_rowdata <- function(x) {
  UseMethod("get_rowdata")
}

#' @rdname get_rowdata
#' @export
get_rowdata.SingleCellExperiment <- function(x) {
  rowData(x) %>% as.data.frame()
}

#' @rdname get_rowdata
#' @export
get_rowdata.seurat <- function(x) {
  data.frame(symbol = rownames(x))
}

#' @rdname get_rowdata
#' @export
get_rowdata.Seurat <- function(x) {
  data.frame(symbol = rownames(x))
}
