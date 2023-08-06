#' @export
get_matrix <- function(x, layer=NULL) {
  if (!is.null(layer))
    X <- x$layers[layer]
  else
    X <- x$X
  X <- Matrix::t(X)
  colnames(X) <- x$obs_names$values
  rownames(X) <- x$var_names$values
　as(X, "CsparseMatrix")
  X
}

#' @export
get_obsm <- function(x, obsm) {
  m <- x$obsm[obsm]
  rownames(m) <- x$obs_names$values

  keys = sub("X_", "", obsm)
  colnames(m) <- paste0(keys, "_", seq_len(ncol(m)))
  m
}

#' @export
create_assay <- function(x, layer=NULL, mod=NULL) {
  if (!is.null(mod))
    x = x[mod]

  X = get_matrix(x, layer=layer)

  SeuratObject::CreateAssayObject(X)
}

#' @export
create_chromatin_assay <- function(x, layer=NULL, mod=NULL, fragments=NULL, annotation=NULL) {
  if (!is.null(mod))
    x = x[mod]

  X = get_matrix(x, layer=layer)
  ranges <- Signac::StringToGRanges(rownames(X), sep = c(":", "-"))
  sel <- as.vector(seqnames(ranges) %in% standardChromosomes(ranges))
  
  Signac::CreateChromatinAssay(X[sel, ], sep = c(":", "-"), fragments=fragments, annotation=annotation)
}

#' @export
create_dimreduc <- function(x, reduction, mod=NULL, assay="RNA") {
  if (!is.null(mod))
    x = x[mod]

  reduction <- get_obsm(x, reduction)
  SeuratObject::CreateDimReducObject(reduction, assay=assay)
}