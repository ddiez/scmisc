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
