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

  lapply(ids, limma::kegga, species = org) |>
    bind_rows(.id = "cluster")
}

#' plot_enrichment
#'
#' Plot the result of run_enrichment function as a heatmap.
#'
#' @param x data.frame with results from run_enrichment.
#' @param group.by grouping variable (default: cluster).
#' @param direction which direction ("Up", "Down") to plot for enrichment (default: c("Up", "Down")).
#' @param FDR cutoff to select significantly enriched terms/pathways.
#' @param col_up colors to use for Up and passed down to ComplexHeatmap::Heatmap.
#' @param col_down colors to use for Down and passed down to ComplexHeatmap::Heatmap.
#' @param top_ann top annotations passed down to ComplexHeatmap::Heatmap.
#' @param circle.scale scale to be applied to the circles.
#' @param ... additional arguments passed down to ComplexHeatmap::Heatmap.
#' @param cluster_columns whether to cluster columns.
#' @param border whether to add border.
#' @param padding padding around heatmap passed down to ComplexHeatmap::draw.
#'
#' @export
plot_enrichment <- function(x, group.by = "cluster", direction = c("Up", "Down"), FDR = 0.01, col_up = c("white", "red"), col_down = c("white", "blue"), top_ann = NULL, circle.scale = .2, cluster_columns = FALSE, border = TRUE, padding = NULL, ...) {
  pathways <- x |> filter(.data[["adj.P.Up"]] < FDR | .data[["adj.P.Down"]] < FDR) |> pull("Pathway")

  direction <- setNames(direction, direction)

  x.color <- lapply(direction, function(n) {
    size.by <- n
    color.by <- paste0("P.", n)
    x |> filter(.data[["Pathway"]] %in% pathways) |>
      select("Pathway", !!group.by, !!color.by) |>
      mutate(color = -log10(.data[[color.by]])) |>
      pivot_wider("Pathway", names_from = group.by, values_from = "color") |>
      column_to_rownames("Pathway") |>
      as.matrix()
  })

  x.size <- lapply(direction, function(n) {
    size.by <- n
    color.by <- paste0("P.", n)
    x |> filter(.data[["Pathway"]] %in% pathways) |>
      select("Pathway", !!group.by, !!size.by) |>
      mutate(size = .data[[size.by]]) |>
      pivot_wider("Pathway", names_from = group.by, values_from = "size") |>
      column_to_rownames("Pathway") |>
      as.matrix()
  })

  m.color <- matrix(NA_real_, ncol = ncol(x.color$Up) * 2, nrow = nrow(x.color$Up))
  colnames(m.color) <- paste0(rep(colnames(x.color$Up), each = 2), "_", rep(c("Up", "Down"), ncol(x.color$Up)))
  rownames(m.color) <- rownames(x.color$Up)
  m.color[, seq(1, ncol(x.color$Up) * 2, by = 2 )] <- x.color$Up
  m.color[, seq(1, ncol(x.color$Up) * 2, by = 2 ) + 1] <- x.color$Down

  m.size <- matrix(NA_real_, ncol = ncol(x.color$Up) * 2, nrow = nrow(x.size$Up))
  colnames(m.size) <- paste0(rep(colnames(x.size$Up), each = 2), "_", rep(c("Up", "Down"), ncol(x.size$Up)))
  rownames(m.size) <- rownames(x.size$Up)
  m.size[, seq(1, ncol(x.size$Up) * 2, by = 2 )] <- x.size$Up
  m.size[, seq(1, ncol(x.size$Up) * 2, by = 2 ) + 1] <- x.size$Down
  m.size <- m.size / max(m.size, na.rm = TRUE)

  split <- sub("_.*", "", colnames(m.color))
  colnames(m.color) <- sub(".*_", "", colnames(m.color))
  col_fun_up <- circlize::colorRamp2(range(m.color), col_up)
  col_fun_down <- circlize::colorRamp2(range(m.color), col_down)
  col_fun <- circlize::colorRamp2(range(m.color), c("white", "black"))

  cell_fun <- function(j, i, x, y, width, height, fill) {
    fill <- ifelse(colnames(m.color)[j] == "Up", col_fun_up(m.color[i, j]), col_fun_down(m.color[i, j]))
    grid::grid.circle(
      x = x,
      y = y,
      r = m.size[i, j] * circle.scale, #* min(unit.c(width, height)),
      gp = grid::gpar(fill = fill, col = NA)
    )
  }

  lgd_up <- ComplexHeatmap::Legend(title = "Up", col_fun = col_fun_up)
  lgd_down <- ComplexHeatmap::Legend(title = "Down", col_fun = col_fun_down)

  h <- ComplexHeatmap::Heatmap(m.color,
                          name = paste0("-log10(P.value)"),
                          cell_fun = cell_fun,
                          col = col_fun,
                          rect_gp = grid::gpar(type = "none"),
                          column_split = split,
                          cluster_columns = cluster_columns,
                          show_heatmap_legend = FALSE,
                          border = border,
                          top_annotation = top_ann, ...)
  if (is.null(padding))
    ComplexHeatmap::draw(h, annotation_legend_list = list(lgd_up, lgd_down), annotation_legend_side = "left")
  else
    ComplexHeatmap::draw(h, annotation_legend_list = list(lgd_up, lgd_down), annotation_legend_side = "left", padding = padding)
}