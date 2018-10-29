#' plot_purity
#'
#' @param x a suitable object.
#' @param col.x name of column to be use in x-axis.
#' @param col.y name of column to be use in y-axis.
#' @param ... parameters passed down to methods.
#'
#' @export
plot_purity <- function(x, ...) {
  UseMethod("plot_purity")
}

#' @rdname plot_purity
#' @export
plot_purity.SingleCellExperiment <- function(x, col.x = "celltype", col.y = "cluster", ...) {
  d <- colData(x) %>% as.data.frame()

  plot_purity(d, col.x = col.x, col.y = col.y, ...)
}

#' @rdname plot_purity
#' @export
plot_purity.data.frame <- function(x, col.x = "celltype", col.y = "cluster", ...) {
  d <- select(x, .data[[col.x]], .data[[col.y]]) %>%
    group_by_(col.y) %>%
    mutate(total = n()) %>%
    group_by_(col.x, col.y) %>%
    summarize(count = n(), percentage = count / first(.data$total))

  ggplot(d, aes_string(col.x, col.y, fill = "percentage")) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red", limits = c(0, 1))
}
