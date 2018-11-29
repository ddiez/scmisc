#' plot_coord
#'
#' @param x object from which to plot coordinates.
#' @param size size of points in geom_point()
#' @param color color or column to color by.
#' @param label name of column used to label plot (e.g. clusters).
#' @param expand logical; whether to expand one column to show presence/absence.
#' @param ... further arguments passed down to get_coord().
#'
#' @export
plot_coord <- function(x, ...) {
  UseMethod("plot_coord")
}

#' @rdname plot_coord
#' @export
plot_coord.seurat <- function(x, size = .1, color = NULL, label = NULL, expand = NULL, ...) {
  d <- get_coord(x, ...)
  plot_coord(d, size = size, color = color, label = label, expand = expand, ...)
}

#' @rdname plot_coord
#' @export
plot_coord.SingleCellExperiment <- function(x, size = .1, color = NULL, label = NULL, expand = NULL, ...) {
  d <- get_coord(x, ...)
  plot_coord(d, size = size, color = color, label = label, expand = expand, ...)
}

#' @rdname plot_coord
#' @export
plot_coord.data.frame <- function(x, size = .1, color = NULL, label = NULL, expand = NULL, ...) {
  d <- x
  if (!is.null(expand)) {
    d <- d %>%
      expand_column(expand) %>%
      arrange_("value")
    p <- ggplot(d, aes_string("dim1", "dim2", color = "value")) +
      geom_point(size = size) +
      scale_color_manual(values = c("grey", "red")) +
      facet_wrap(~.data[[expand]])
  } else {
    if (!is.null(color))
      d <- d %>% arrange_(color)

    p <- ggplot(d, aes_string("dim1", "dim2")) +
      geom_point(size = size)

    if (!is.null(color)) {
      p <- p + aes_string(color = color)
      if (is.numeric(d[[color]]))
        p <- p + scale_color_gradient(low = "grey", high = "red")
    }

    if (!is.null(label)) {
      dd <- d %>% group_by_(label) %>%
        summarize(dim1 = mean(.data$dim1), dim2 = mean(.data$dim2))
      p <- p + geom_text(aes_string(label = label), data = dd, color = "black")
    }
  }

  p
}

#' plot_purity
#'
#' Plot heatmap of purity index of specified columns.
#'
#' @param x a suitable object.
#' @param col.x name of column to be use in x-axis.
#' @param col.y name of column to be use in y-axis.
#' @param label logical; whether to add rounded values of purity to tiles.
#' @param label.size size of labels text.
#' @param ... parameters passed down to methods.
#'
#' @export
plot_purity <- function(x, ...) {
  UseMethod("plot_purity")
}

#' @rdname plot_purity
#' @export
plot_purity.SingleCellExperiment <- function(x, ...) {
  plot_purity(colData(x), ...)
}

#' @rdname plot_purity
#' @export
plot_purity.DataFrame <- function(x, ...) {
  plot_purity(as.data.frame(x), ...)
}

#' @rdname plot_purity
#' @export
plot_purity.data.frame <- function(x, col.x, col.y, label = FALSE, label.size = 5, ...) {
  if (missing(col.x)) stop("specify name of column for x-axis.")
  if (missing(col.y)) stop("specify name of column for y-axis.")

  d <- compute_purity(x, col.x, col.y)

  p <- ggplot(d, aes_string(col.x, col.y, fill = "purity")) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red", limits = c(0, 1))

  if (label) {
    d <- d %>% mutate(purity = format(round(.data$purity, 2)))
    p <- p + geom_text(aes_string(col.x, col.y, label = "purity"), data = d, size = label.size, inherit.aes = FALSE)
  }

  p
}

#' plot_jaccard
#'
#' Plot heatmap of jaccard index of specified columns.
#'
#' @param x a suitable object.
#' @param col.x name of column to be use in x-axis.
#' @param col.y name of column to be use in y-axis.
#' @param label logical; whether to add rounded values of jaccard index to tiles.
#' @param label.size size of labels text.
#' @param ... parameters passed down to methods.
#'
#' @export
plot_jaccard <- function(x, ...) {
  UseMethod("plot_jaccard")
}

#' @rdname plot_jaccard
#' @export
plot_jaccard.SingleCellExperiment <- function(x, ...) {
  plot_jaccard(colData(x), ...)
}

#' @rdname plot_jaccard
#' @export
plot_jaccard.DataFrame <- function(x, ...) {
  plot_jaccard(as.data.frame(x), ...)
}

#' @rdname plot_jaccard
#' @export
plot_jaccard.data.frame <- function(x, col.x, col.y, label = FALSE, label.size = 5, ...) {
  if (missing(col.x)) stop("specify name of column for x-axis.")
  if (missing(col.y)) stop("specify name of column for y-axis.")

  d <- compute_jaccard(x, col.x, col.y)

  p <- ggplot(d, aes_string(col.x, col.y, fill = "jaccard")) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red", limit = c(0, 1))

  if (label) {
    d <- d %>% mutate(jaccard = format(round(.data$jaccard, 2)))
    p <- p + geom_text(aes_string(col.x, col.y, label = "jaccard"), data = d, size = label.size, inherit.aes = FALSE)
  }

  p
}
