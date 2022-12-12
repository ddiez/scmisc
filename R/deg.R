#' top_deg
#'
#' @param x data.frame or similar with DEGs.
#' @param fdr FDR used for filtering.
#' @param lfc log2 fold change used for filtering.
#' @param n number of genes per group to return.
#' @param group column name to group results.
#'
#' @export
top_deg <- function(x, fdr = 0.01, lfc = 1, n = 10, group = "cluster") {
  top <- x |> filter(abs(.data[["avg_log2FC"]]) > lfc, .data[["p_val_adj"]] < fdr) |> arrange(.data[["p_val_adj"]])
  top |> group_by(.data[[group]]) |> slice_head(n = n)
}

#' plot_deg_barplot
#'
#' Plot top DEG.
#'
#' @param x data.frame or similar with DEGs.
#' @param fdr FDR used for filtering.
#' @param lfc log2 fold change used for filtering.
#' @param group column name to group results.
#'
#' @export
plot_deg_barplot <- function(x, lfc=1, fdr=0.01, group="cluster") {
  d <- x |> filter(p_val_adj < fdr, abs(avg_log2FC) > lfc) |>
    mutate(direction=ifelse(avg_log2FC>0, "Up", "Down")) |>
    mutate(direction=factor(.data[["direction"]], levels=c("Up", "Down"))) |>
    count(.data[[group]], .data[["direction"]])
  ggplot(d, aes(.data[[group]], .data[["n"]], fill=.data[["direction"]])) +
    geom_col(position="dodge") +
    scale_fill_manual(values=list("Up"="red", "Down"="blue")) +
    Seurat::RotatedAxis() +
    labs(x="", y="")
}
