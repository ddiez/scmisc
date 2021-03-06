% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coord.R
\name{get_coord}
\alias{get_coord}
\alias{get_coord.SingleCellExperiment}
\alias{get_coord.CellDataSet}
\alias{get_coord.seurat}
\alias{get_coord.Seurat}
\title{get_coord}
\usage{
get_coord(x, coord.name = NULL, add.cols = TRUE, add.exprs = NULL, ...)

\method{get_coord}{SingleCellExperiment}(x, coord.name = NULL, add.cols = TRUE, add.exprs = FALSE, ...)

\method{get_coord}{CellDataSet}(x, coord.name = "A", add.cols = TRUE, add.exprs = FALSE, ...)

\method{get_coord}{seurat}(x, coord.name = "tsne", add.cols = TRUE, add.exprs = FALSE, ...)

\method{get_coord}{Seurat}(
  x,
  coord.name = NULL,
  add.cols = TRUE,
  add.exprs = FALSE,
  assay = NULL,
  slot = "data",
  ...
)
}
\arguments{
\item{x}{object to obtain coordinates.}

\item{coord.name}{name of the holder for coordinates.}

\item{add.cols}{logical; whether to annotate each coordinate with column data (i.e. colData).}

\item{add.exprs}{if not NULL a list of genes to add counts obtained with get_expression().}

\item{...}{arguments passed down to methods.}

\item{assay}{name of assay for Seurat objects.}

\item{slot}{name of slot for Seurat objects.}
}
\value{
A data.frame object.
}
\description{
Extract coordinates (e.g. from dimensionality reduction) from bio datasets,
optionally with annotations to help plotting them with ggplot2 package.
}
