#' run_enrichment
#'
#' @param x an object with DEGs.
#' @param type which ontology to use for enrichment (default: KEGG).
#' @param group.by which variable to use for grouping (default: cluster).
#' @param use.column which column contains entrez gene identifiers (default: entrezgene).
#' @param org organisms (default: Mm).
#' @param FDR cutoff for selection of DEG used for enrichment.
#' @param p.adjust.method method to adjust p.values from enrichment (default: bonferroni).
#' @param ... arguments passed down to methods.
#'
#' @export
run_enrichment <- function(x, ...) {
  UseMethod("run_enrichment")
}

#' @rdname run_enrichment
#' @export
run_enrichment.data.frame <- function(x, type = "kegg", group.by = "cluster", use.column = "entrezgene", org = "Mm", FDR = 0.01, p.adjust.method = "bonferroni", ...) {
  if (!is.null(FDR))
    x <- x |> filter(.data[["p_val_adj"]] <= FDR)

  res <- kegg_enrichment(x, group.by = group.by, use.column = use.column, org = org)

  res |> mutate(
    adj.P.Up = p.adjust(.data[["P.Up"]], method = p.adjust.method),
    adj.P.Down = p.adjust(.data[["P.Down"]], method = p.adjust.method),
  )
}

kegg_enrichment <- function(x, group.by = "cluster", use.column = "entrezgene", org = "Mm") {
  x <- x |> drop_na(any_of(use.column))

  groups <- x[[group.by]]
  if (is.factor(groups))
    groups <- levels(groups)
  else
    groups <- unique(groups)

  ids <- lapply(groups, function(group) {
    genes <- list(
      Up = x |> filter(.data[[group.by]] == !!group, .data[["avg_log2FC"]] > 0) |> pull(use.column),
      Down = x |> filter(.data[[group.by]] == !!group, .data[["avg_log2FC"]] < 0) |> pull(use.column))
  })
  names(ids) <- groups

  lapply(ids, limma::kegga, species = org) |> bind_rows(.id = "cluster")
}

#' plot_enrichment
#'
#' Plot the result of run_enrichment function as a heatmap.
#'
#' @param x data.frame with results from run_enrichment.
#' @param group.by grouping variable (default: cluster).
#' @param color.by color variable (default: P.Up).
#' @param FDR cutoff to select significantly enriched terms/pathways.
#' @param col colors to use and passed down to ComplexHeatmap::Heatmap.
#' @param top_ann top annotations passed down to ComplexHeatmap::Heatmap.
#' @param ... additional arguments passed down to ComplexHeatmap::Heatmap.
#'
#' @export
plot_enrichment <- function(x, group.by = "cluster", color.by = "P.Up", FDR = 0.01, col = c("white", "blue"), top_ann = NULL, ...) {
  pathways <- x |> filter(.data[["adj.P.Up"]] < FDR | .data[["adj.P.Down"]] < FDR) |> pull("Pathway")

  m <- x |> filter(.data[["Pathway"]] %in% pathways) |>
    select("Pathway", !!group.by, !!color.by) |>
    mutate(color = -log10(.data[[color.by]])) |>
    pivot_wider("Pathway", names_from = group.by, values_from = "color") |>
    column_to_rownames("Pathway") |>
    as.matrix()

  ComplexHeatmap::Heatmap(m, col = col, top_annotation = top_ann, ...)
}
