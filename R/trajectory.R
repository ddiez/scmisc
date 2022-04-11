
#' plot_trajectory_graph
#'
#' Plot the MST graph used to estimate trajectories.
#'
#' @param x and object with trajectory information.
#' @param layout layout for the graph.
#' @param ... arguments passed down to methods.
#'
#' @export
plot_trajectory_graph <- function(x, ...) {
  UseMethod("plot_trajectory_graph")
}

#' @rdname plot_trajectory_graph
#' @export
plot_trajectory_graph.SingleCellExperiment <- function(x, ...) {
  plot_trajectory_graph(slingshot::SlingshotDataSet(x), ...)
}

#' @rdname plot_trajectory_graph
#' @export
plot_trajectory_graph.SlingshotDataSet <- function(x, ...) {
  g <- slingshot::slingMST(x)
  g <- igraph::graph_from_adjacency_matrix(g, mode="undirected")
  g <- tidygraph::as_tbl_graph(g, directed = FALSE)
  g <- g |> tidygraph::activate("nodes") |>
    mutate(cluster = "middle")

  sc <- slingshot::slingParams(x)[["start.clus"]]
  if (!is.null(sc)) {
    g <- g |> tidygraph::activate("nodes") |>
      mutate(cluster = case_when(
        .data$name %in% sc ~ "start",
        TRUE ~ .data$cluster
      ))
  }

  ec <- slingshot::slingParams(x)[["end.clus"]]
  if (!is.null(ec)) {
    g <- g |> tidygraph::activate("nodes") |>
      mutate(cluster = case_when(
        name %in% ec ~ "end",
        TRUE ~ .data$cluster
      ))
  }

  plot_trajectory_graph(g, ...)
}

#' @rdname plot_trajectory_graph
#' @export
plot_trajectory_graph.tbl_graph <- function(x, layout = "nicely", ...) {
  l <- ggraph::create_layout(x, layout = layout)
  plot_trajectory_graph(l)
}

#' @rdname plot_trajectory_graph
#' @export
plot_trajectory_graph.layout_ggraph <- function(x, ...) {
  ggraph::ggraph(x) +
    ggraph::geom_node_text(aes_string(label = "name", color = "cluster")) +
    ggraph::geom_edge_fan(
      end_cap = ggraph::circle(10, "points"),
      start_cap = ggraph::circle(10, "points"),
    ) +
    ggraph::theme_graph()
}

#' plot_trajectory
#'
#' Plot the cell embedings with colors indicating the increasing trajectory pseudotime.
#'
#' @param x  and object with trajectory information.
#' @param coord.name name of the reduced dimension to use.
#' @param ... arguments passed down to methods.
#'
#' @export
plot_trajectory <- function(x, ...) {
  UseMethod("plot_trajectory")
}

#' @rdname plot_trajectory
#' @export
plot_trajectory.SingleCellExperiment <- function(x, coord.name = NULL, ...) {
  if (is.null(coord.name))
    coord.name = reducedDimNames(x)[1]

  d <- get_coord(x, coord.name) |>
    gather("curve", "pseudotime", starts_with("slingPseudotime"))

  ggplot(d, aes_string("dim1", "dim2", color = "pseudotime")) +
    geom_point(size = .5) +
    scale_color_distiller(palette = "Spectral") +
    facet_wrap(~curve)
}

#' @rdname plot_trajectory
#' @export
plot_trajectory.SlingshotDataSet <- function(x, ...) {
  d <- reducedDim(x) |> as_tibble() |>
    rename(dim1 = 1, dim2 = 2)

  d <- cbind(d, slingshot::slingPseudotime(x)) |>
    gather("curve", "pseudotime", starts_with("curve"))

  ggplot(d, aes_string("dim1", "dim2", color = "pseudotime")) +
    geom_point(size = .5) +
    scale_color_distiller(palette = "Spectral") +
    facet_wrap(~curve)
}

#' plot_pseudotime_gene
#'
#'
#' @param x  and object with pseudotime information.
#' @param ... arguments passed down to methods.
#'
#' @export
plot_pseudotime_gene <- function(x, ...) {
  UseMethod("plot_pseudotime_gene")
}

#' @rdname plot_pseudotime_gene
#' @export
plot_pseudotime_gene.cell_data_set <- function(x, features=NULL, cutoff=0, combine=TRUE, assay="logcounts", ...) {
  features <- intersect(features, rownames(x))
  exprs <- SummarizedExperiment::assay(x, assay)[features, , drop=FALSE]
  exprs <- Matrix::t(exprs)
  d <- cbind(get_coord(x), exprs)
  d$pseudotime <- monocle3::pseudotime(x)

  p <- lapply(features, function(feature) {
    ggplot(d |> filter(.data[[feature]] > cutoff), aes(.data[["pseudotime"]], .data[[feature]], color = .data[["pseudotime"]])) +
      geom_jitter(width=.5, size=.1) +
      scale_color_viridis_c(option="plasma") +
      geom_smooth(method="lm", formula=y~ splines::ns(x, df=3), se=FALSE, color="violetred")
  })

  if(length(p) == 1) return(p[[1]])

  if (combine)
    p |> wrap_plots()
  else
    p
}

#' plot_pseudotime_heatmap
#'
#'
#' @param x  and object with pseudotime information.
#' @param ... arguments passed down to methods.
#'
#' @export
plot_pseudotime_heatmap <- function(x, ...) {
  UseMethod("plot_pseudotime_heatmap")
}

#' @rdname plot_pseudotime_heatmap
#' @export
plot_pseudotime_heatmap.Seurat <- function(x, features, assay="RNA", slot="data", reduction="pseudotime", ...) {
  pseudotime <- Embeddings(x, reduction=reduction)[, 1]

  sel.good <- ! is.infinite(pseudotime)
  pseudotime <- pseudotime[sel.good]
  x <- x[, sel.good]

  pseudotime <- sort(pseudotime)
  cells <- names(pseudotime)

  m <- GetAssayData(x, assay=assay, slot=slot)
  m <- m[features, cells]|> as.matrix()
  m <- t(scale(t(m)))

  #pseudo_color <- circlize::colorRamp2(range(pseudotime, na.rm=TRUE), c("white", "red"))
  r <- range(pseudotime, na.rm=TRUE)
  pseudo_color <- circlize::colorRamp2(seq(r[1], r[2], length.out=5), viridis::magma(5))

  #r <- range(m, na.rm=TRUE)
  #m_color <- circlize::colorRamp2(seq(r[1], r[2], length.out=5), viridis::magma(5))

  df <- data.frame(pseudotime = pseudotime[cells])
  top_cols <- list(pseudotime=pseudo_color)
  top_ann <- ComplexHeatmap::columnAnnotation(df=df, col=top_cols)

  plot_heatmap(m, cluster_columns=FALSE, top_ann=top_ann, ...)
}

#' plot_pseudotime_modules
#'
#'
#' @param x  and object with pseudotime information.
#' @param ... arguments passed down to methods.
#'
#' @export
plot_pseudotime_modules <- function(x, ...) {
  UseMethod("plot_pseudotime_modules")
}

#' @rdname plot_pseudotime_modules
#' @export
plot_pseudotime_modules.Seurat <- function(x, gene_modules, reduction="pseudotime", assay="RNA", slot="data", filter.zero=TRUE, add.jitter=FALSE, width=1) {
  pseudotime <- Embeddings(x, reduction=reduction)[, 1]
  pseudotime <- pseudotime[!is.infinite(pseudotime)]


  pseudotime <- sort(pseudotime)
  cells <- names(pseudotime)

  m <- GetAssayData(x, assay=assay, slot=slot)
  m <- m[, cells]

  modules <- unique(gene_modules$module)
  d <- lapply(modules, function(module) {
    genes <- gene_modules |> filter(.data[["module"]] == !!module) |> pull(gene)
    means <- colMeans(m[genes, ])
    data.frame(index=seq_along(pseudotime), pseudotime=pseudotime, mean=means, module=module)
  }) |> bind_rows()

  if (filter.zero)
    d <- d |> filter(mean>0)

  if (add.jitter)
    points <- geom_jitter(size = .1, width=width)
  else
    points <- geom_point(size=.1)

  ggplot(d, aes(pseudotime, mean, color=pseudotime)) +
    points +
    geom_smooth(color="violetred", method="gam", formula=y ~ s(x, bs = "cs")) +
    scale_color_viridis_c(option="magma") +
    facet_wrap(~module, scales="free_y")
}
