#' plot_coord
#'
#' @param x object from which to plot coordinates.
#' @param size size of points in geom_point()
#' @param color color or column to map to color..
#' @param shape shape or column to map to shape.
#' @param label name of column used to label plot (e.g. clusters).
#' @param expand logical; whether to expand one column to show presence/absence.
#' @param ncol number of columns used in facet_wrap.
#' @param ... further arguments passed down to get_coord().
#'
#' @export
plot_coord <- function(x, ...) {
  UseMethod("plot_coord")
}

#' @rdname plot_coord
#' @export
plot_coord.seurat <- function(x, size = 1, color = NULL, shape = NULL, label = NULL, expand = NULL, ncol = NULL, ...) {
  d <- get_coord(x, ...)
  plot_coord(d, size = size, color = color, shape = shape, label = label, expand = expand, ncol = ncol, ...)
}

#' @rdname plot_coord
#' @export
plot_coord.Seurat <- function(x, size = 1, color = NULL, shape = NULL, label = NULL, expand = NULL, ncol = NULL, ...) {
  d <- get_coord(x, ...)
  plot_coord(d, size = size, color = color, shape = shape, label = label, expand = expand, ncol = ncol, ...)
}

#' @rdname plot_coord
#' @export
plot_coord.SingleCellExperiment <- function(x, size = 1, color = NULL, shape = NULL, label = NULL, expand = NULL, ncol = NULL, ...) {
  d <- get_coord(x, ...)
  plot_coord(d, size = size, color = color, shape = shape, label = label, expand = expand, ncol = ncol, ...)
}

#' @rdname plot_coord
#' @export
plot_coord.data.frame <- function(x, size = 1, color = NULL, shape = NULL, label = NULL, expand = NULL, ncol = NULL, ...) {
  d <- x
  if (!is.null(expand)) {
    d <- d |>
      expand_column(expand) |>
      arrange(.data[["value"]])

    if (length(expand) == 1) {
      p <- ggplot(d, aes_string("dim1", "dim2", color = "value")) +
        geom_point(size = size) +
        scale_color_manual(values = c("grey", "red")) +
        facet_wrap(~.data[[expand]], ncol = ncol) +
        guides(color = "none")
    }

    if (length(expand) == 2) {
      p <- ggplot(d, aes_string("dim1", "dim2", color = "value")) +
        geom_point(size = size) +
        scale_color_manual(values = c("grey", "red")) +
        facet_grid(rows = vars(.data[[expand[2]]]), cols = vars(.data[[expand[1]]])) +
        guides(color = "none")
    }
  } else {
    if (!is.null(color))
      d <- d |> arrange(.data[[color]])

    p <- ggplot(d, aes_string("dim1", "dim2")) +
      geom_point(size = size)

    if (!is.null(color)) {
      p <- p + aes_string(color = color)
      if (is.numeric(d[[color]]))
        p <- p + scale_color_gradient(low = "grey", high = "red")
    }

    if (!is.null(shape)) {
      p <- p + aes_string(shape = shape)
    }

    if (!is.null(label)) {
      dd <- d |> group_by_(label) |>
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
#' @param drop whether to drop unused axis levels.
#' @param ... parameters passed down to methods.
#'
#' @export
plot_purity <- function(x, ...) {
  UseMethod("plot_purity")
}

#' @rdname plot_purity
#' @export
plot_purity.Seurat <- function(x, col.x, col.y, ...) {
  plot_purity(FetchData(x, vars = c(col.x, col.y)), col.x = col.x, col.y = col.y, ...)
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
plot_purity.data.frame <- function(x, col.x, col.y, label = FALSE, label.size = 5, drop = FALSE, ...) {
  if (missing(col.x)) stop("specify name of column for x-axis.")
  if (missing(col.y)) stop("specify name of column for y-axis.")

  d <- compute_purity(x, col.x, col.y)

  p <- ggplot(d, aes_string(col.x, col.y, fill = "purity")) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red", limits = c(0, 1))

  p <- p + scale_x_discrete(drop = FALSE, expand = c(0, 0)) +
    scale_y_discrete(drop = FALSE, expand = c(0, 0))

  if (label) {
    d <- d |> mutate(purity = format(round(.data$purity, 2)))
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
#' @param drop whether to drop unused axis levels.
#' @param ... parameters passed down to methods.
#'
#' @export
plot_jaccard <- function(x, ...) {
  UseMethod("plot_jaccard")
}

#' @rdname plot_jaccard
#' @export
plot_jaccard.Seurat <- function(x, col.x, col.y, ...) {
  plot_jaccard(FetchData(x, vars = c(col.x, col.y)), col.x = col.x, col.y = col.y, ...)
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
plot_jaccard.data.frame <- function(x, col.x, col.y, label = FALSE, label.size = 5, drop = FALSE, ...) {
  if (missing(col.x)) stop("specify name of column for x-axis.")
  if (missing(col.y)) stop("specify name of column for y-axis.")

  d <- compute_jaccard(x, col.x, col.y)

  p <- ggplot(d, aes_string(col.x, col.y, fill = "jaccard")) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red", limit = c(0, 1))

  p <- p + scale_x_discrete(drop = FALSE, expand = c(0, 0)) +
      scale_y_discrete(drop = FALSE, expand = c(0, 0))

  if (label) {
    d <- d |> mutate(jaccard = format(round(.data$jaccard, 2)))
    p <- p + geom_text(aes_string(col.x, col.y, label = "jaccard"), data = d, size = label.size, inherit.aes = FALSE)
  }

  p
}

#' plot_violin
#'
#' Makes violin plot, by default using ggbeeswarm package.
#'
#' @param x object to plot data from.
#' @param feature name of the feature to plot.
#' @param group grouping variable.
#' @param size size of the points.
#' @param ... further arguments passed down to get_expression().
#'
#' @export
plot_violin <- function(x, feature, group, ...) {
  UseMethod("plot_violin")
}

#' @rdname plot_violin
#' @export
plot_violin.Seurat <- function(x, feature, group, size = .1, ...) {
  d <- cbind(x@meta.data, expression = get_expression(x, feature, ...))
  plot_violin(d, feature = feature, group = group, size = size)
}

#' @rdname plot_violin
#' @export
plot_violin.SingleCellExperiment <- function(x, feature, group, size = .1, ...) {
  d <- cbind(as.data.frame(SummarizedExperiment::colData(x)), expression = get_expression(x, feature, ...))
  plot_violin(d, feature = feature, group = group, size = size)
}

#' @rdname plot_violin
#' @export
plot_violin.data.frame <- function(x, feature, group, size = .1, ...) {
  ggplot(x, aes_string(group, "expression")) + ggbeeswarm::geom_quasirandom(groupOnX = TRUE, size = size) +
    labs(title = feature)
}


#' plot_heatmap
#'
#' @param x object to plot.
#' @param assay assay name for Seurat and SingleCellExperiment objects.
#' @param slot slot name for Seurat objects.
#' @param scale whether to scale data data.
#' @param top_ann names of columns to be used as top annotations.
#' @param top_ann_col color definition for the categories in the top annotations.
#' @param show_column_names whether to show column names (default: FALSE).
#' @param ... arguments passed down to ComplexHeatmap::Heatmap.
#'
#' @export
plot_heatmap <- function(x, ...) {
  UseMethod("plot_heatmap")
}

#' @rdname plot_heatmap
#' @export
plot_heatmap.Seurat <- function(x, assay = NULL, slot = "data", scale = TRUE, top_ann = NULL, top_ann_col = NULL, ...) {
  if (!is.null(top_ann)) {
    df <- x[[]][, top_ann, drop = FALSE]
    top_ann <- ComplexHeatmap::columnAnnotation(df = df, col = top_ann_col)
  }

  if (slot == "scale.data")
    scale = FALSE

  x <- GetAssayData(x, assay = assay, slot = slot)
  plot_heatmap(as.matrix(x), scale = scale, top_ann = top_ann, ...)
}


#' @rdname plot_heatmap
#' @export
plot_heatmap.SingleCellExperiment <- function(x, assay = "logcounts", top_ann = NULL, top_ann_col = NULL, ...) {
  if (!is.null(top_ann)) {
    df <- SummarizedExperiment::colData(x)[, top_ann]
    top_ann <- ComplexHeatmap::columnAnnotation(df = df, col = top_ann_col)
  }

  x <- SummarizedExperiment::assay(x, assay)
  plot_heatmap(as.matrix(x), top_ann = top_ann, ...)
}


#' @rdname plot_heatmap
#' @export
plot_heatmap.matrix <- function(x, scale = TRUE, show_column_names = FALSE, ...) {
  if (scale)
    x <- t(scale(t(x)))

  ComplexHeatmap::Heatmap(x, name = "expression", show_column_names = show_column_names, ...)
}


#' plot_loadings
#'
#' @param x object to plot.
#' @param features name of features to plot.
#' @param reduction name of reduction to plot.
#' @param cluster_columns whether to cluster heatmap columns.
#' @param ... arguments passed down to ComplexHeatmap::Heatmap.
#'
#' @export
plot_loadings <- function(x, features, reduction = "pca", cluster_columns = FALSE, ...) {
  UseMethod("plot_loadings")
}

#' @rdname plot_loadings
#' @export
plot_loadings.Seurat <- function(x, features, reduction = "pca", cluster_columns = FALSE, ...) {
  d <- SeuratObject::Loadings(x, reduction = reduction)
  features <- intersect(features, rownames(d))
  ComplexHeatmap::Heatmap(d[features, ], cluster_columns = cluster_columns, name = "loadings", ...)
}

#' plot_velocity
#'
#' @param x data.frame with velocity information.
#' @param meta cell information.
#' @param color column from meta to color cells.
#' @param grid whether to plot vector field.
#' @param arrow.color color for arrows.
#' @param arrow.size size for arrows.
#'
#' @export
plot_velocity <- function(x, meta = NULL, color = NULL, grid = FALSE, arrow.color = "black", arrow.size = .2) {
  d <- x$arrows |>
    as.data.frame() |>
    rownames_to_column("barcode")

  if (! is.null(meta)) {
    d <- d |> left_join(meta |> rownames_to_column("barcode"), by = "barcode")
  }

  if (grid) {
    da <- x$garrow |>
      as.data.frame()
  } else {
    da <- d
  }

  if (!is.null(meta) && !is.null(color)) {
    p <- ggplot(d, aes(x = .data[["x0"]], y = .data[["y0"]], color = .data[[color]])) +
      geom_point()
  } else {
    p <- ggplot(d, aes(x = .data[["x0"]], y = .data[["y0"]])) +
      geom_point(color = "grey")
  }
  p + geom_segment(aes(x = .data[["x0"]], y = .data[["y0"]], xend = .data[["x1"]], yend = .data[["y1"]]), data = da, arrow = arrow(length = unit(1, "mm"), type = "closed"), color = arrow.color, size = arrow.size) +
    theme_void()
}
