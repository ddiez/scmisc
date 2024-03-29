% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrichment.R
\name{plot_enrichment}
\alias{plot_enrichment}
\title{plot_enrichment}
\usage{
plot_enrichment(
  x,
  group.by = "cluster",
  direction = c("Up", "Down"),
  FDR = 0.01,
  col_up = c("white", "red"),
  col_down = c("white", "blue"),
  top_ann = NULL,
  circle.scale = 0.2,
  cluster_columns = FALSE,
  border = TRUE,
  padding = NULL,
  ...
)
}
\arguments{
\item{x}{data.frame with results from run_enrichment.}

\item{group.by}{grouping variable (default: cluster).}

\item{direction}{which direction ("Up", "Down") to plot for enrichment (default: c("Up", "Down")).}

\item{FDR}{cutoff to select significantly enriched terms/pathways.}

\item{col_up}{colors to use for Up and passed down to ComplexHeatmap::Heatmap.}

\item{col_down}{colors to use for Down and passed down to ComplexHeatmap::Heatmap.}

\item{top_ann}{top annotations passed down to ComplexHeatmap::Heatmap.}

\item{circle.scale}{scale to be applied to the circles.}

\item{cluster_columns}{whether to cluster columns.}

\item{border}{whether to add border.}

\item{padding}{padding around heatmap passed down to ComplexHeatmap::draw.}

\item{...}{additional arguments passed down to ComplexHeatmap::Heatmap.}
}
\description{
Plot the result of run_enrichment function as a heatmap.
}
