
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
  g <- tidygraph::as_tbl_graph(g)
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
