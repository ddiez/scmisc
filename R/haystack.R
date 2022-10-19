#' plot_result_rank
#'
#' @param x haystack object or summary data.frame.
#' @param highlight gene names to highlight.
#' @param sort.by rank the genes by the following column from the summary data.frame.
#'
#' @export
#'
plot_result_rank <- function(x, highlight=NULL, sort.by="log.p.adj") {
  UseMethod("plot_result_rank")
}

#' @rdname plot_result_rank
#' @export
plot_result_rank.haystack <- function(x, highlight=NULL, sort.by="log.p.adj") {
  sum <- show_result_haystack(x)
  plot_result_rank(sum, highlight=highlight, sort.by=sort.by)
}

#' @rdname plot_result_rank
#' @export
plot_result_rank.data.frame <- function(x, highlight=NULL, sort.by="log.p.adj") {
  x <- x |> rownames_to_column("gene")
  x <- x |> arrange(sort.by)
  x <- x |> mutate(index=seq_len(nrow(x)))
  p <- ggplot(x, aes(.data[["index"]], .data[[sort.by]])) +
    geom_point()

  if (!is.null(highlight)) {
    x <- x |> filter(.data[["gene"]] %in% highlight)
    p <- p + ggrepel::geom_text_repel(aes(label=.data[["gene"]]), color="violetred", data=x, min.segment.length=0)
  }
  p
}
