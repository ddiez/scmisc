% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deg.R
\name{plot_deg_barplot}
\alias{plot_deg_barplot}
\title{plot_deg_barplot}
\usage{
plot_deg_barplot(x, lfc = 1, fdr = 0.01, group = "cluster")
}
\arguments{
\item{x}{data.frame or similar with DEGs.}

\item{lfc}{log2 fold change used for filtering.}

\item{fdr}{FDR used for filtering.}

\item{group}{column name to group results.}
}
\description{
Plot top DEG.
}
