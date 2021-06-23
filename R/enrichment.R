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
