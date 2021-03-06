% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter.R
\name{filter_genes}
\alias{filter_genes}
\alias{filter_genes.SingleCellExperiment}
\alias{filter_genes.matrix}
\title{filter_genes}
\usage{
filter_genes(x, ...)

\method{filter_genes}{SingleCellExperiment}(x, min.cells = 10, min.count = 0, assay.name = "counts", ...)

\method{filter_genes}{matrix}(x, min.cells = 10, min.count = 0, ...)
}
\arguments{
\item{x}{object to filter.}

\item{...}{arguments passed down to methods.}

\item{min.cells}{minimum number of cells with min.count for gene to be preserved.}

\item{min.count}{minimum number of counts to be considered expressed.}

\item{assay.name}{name of assay to obtain matrix of counts.}
}
\description{
filter_genes
}
