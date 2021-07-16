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
  colData(x) |> as.data.frame()
}

#' @rdname get_coldata
#' @export
get_coldata.seurat <- function(x) {
  x@meta.data |> as.data.frame()
}

#' @rdname get_coldata
#' @export
get_coldata.Seurat <- function(x) {
  x@meta.data |> as.data.frame()
}


#' get_rowdata
#'
#' @param x object.
#' @param assay name of assay.
#' @param ... arguments passed down to methods.
#'
#' @export
get_rowdata <- function(x, ...) {
  UseMethod("get_rowdata")
}

#' @rdname get_rowdata
#' @export
get_rowdata.SingleCellExperiment <- function(x, ...) {
  rowData(x) |> as.data.frame()
}

#' @rdname get_rowdata
#' @export
get_rowdata.seurat <- function(x, ...) {
  data.frame(symbol = rownames(x))
}

#' @rdname get_rowdata
#' @export
get_rowdata.Seurat <- function(x, assay = NULL, ...) {
  if (is.null(assay)) assay <- DefaultAssay(x)
  data.frame(symbol = rownames(GetAssayData(x, assay = assay)))
}

#' get_data
#'
#' Get miscellaneous data from objects.
#'
#' @param x  object to extract data from.
#' @param coord.name name of coordinates.
#' @param feature.name name of features.
#' @param meta.data whether to include meta.data.
#' @param assay name of assay to use for Seurat objects.
#' @param slot name of slot to use for Seurat objects.
#' @param ... argument passed down to methods.
#'
#' @export
get_data <- function(x, coord.name = NULL, feature.name = NULL, meta.data = TRUE, ...) {
  UseMethod("get_data")
}

#' @rdname get_data
#' @export
get_data.Seurat <- function(x, coord.name = NULL, feature.name = NULL, meta.data = TRUE, assay = NULL, slot = "data", ...) {
  meta <- NULL
  if (meta.data) {
    meta <- x[[]]
  }

  coord <- NULL
  if (!is.null(coord.name)) {
    if (isTRUE(coord.name)) {
      coord.name <- SeuratObject::Reductions(x)[1]
    }
    coord <- get_coord(x, coord.name = coord.name)
  }

  feature <- NULL
  if(!is.null(feature.name)) {
    feature <- sapply(feature.name, function(gene) {
      get_expression(x, gene, assay.name = assay, slot = slot)
    })
  }

  list(embedings = coord, meta = meta, features = feature)
}

