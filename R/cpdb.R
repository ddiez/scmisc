#' read_cpdb_out
#'
#' Read cellphonedb results.
#'
#' @param path path to the folder containing cellphonedb results.
#'
#' @export
#'
read_cpdb_out <- function(path) {
  cols <- c("id_cp_interaction", "interacting_pair", "partner_a", "partner_b", "gene_a", "gene_b", "secreted", "receptor_a", "receptor_b", "annotation_strategy", "is_integrin")
  d_pval <- read_tsv(file.path(path, "pvalues.txt"), show_col_types=FALSE)
  d_mean <- read_tsv(file.path(path, "means.txt"), show_col_types=FALSE)

  cells <- setdiff(colnames(d_pval), cols)
  proteins <- d_pval[["interacting_pair"]]

  m_pval <- d_pval[, cells]
  m_mean <- d_mean[, cells]
  meta <- d_pval[, cols]

  list(mean=d_mean, pval=d_pval, cols=cols, cells=cells, proteins=proteins, meta=meta, m_mean=m_mean, m_pval=m_pval)
}

#' plot_cpdb_dotplot
#'
#' Plot cellphonedb results as a dotplot.
#'
#' @param path path to the folder containing cellphonedb results.
#' @param cutoff cutoff for pvalues.
#' @param filter filter for cell pairs using regular expression.
#' @param cells subset cell pairs.
#' @param proteins subset protein pairs.
#' @param cluster_cells logical; whether to cluster cells.
#' @param cluster_proteins logigcal; whether to cluster proteins.
#'
#' @export
#'
plot_cpdb_dotplot <- function(x, cutoff=1e-3, filter=NULL, cells=NULL, proteins=NULL, cluster_cells=FALSE, cluster_proteins=FALSE) {
  d.pval <- x$pval
  d.mean <- x$mean

  cols <- x$cols

  if (is.null(proteins))
    proteins <- x$proteins

  if (is.null(cells))
    cells <- x$cells

  if (!is.null(filter))
    cells <- grep(filter, cells, value=TRUE)

  d.pval <- d.pval[, c("interacting_pair", cells)] |> filter(.data[["interacting_pair"]] %in% proteins)
  d.mean <- d.mean[, c("interacting_pair", cells)] |> filter(.data[["interacting_pair"]] %in% proteins)

  d.pval[["interacting_pair"]] <- make.names(d.pval[["interacting_pair"]], unique=TRUE)
  d.mean[["interacting_pair"]] <- make.names(d.mean[["interacting_pair"]], unique=TRUE)

  pval <- d.pval |> column_to_rownames("interacting_pair") |> as.matrix()
  sig_inter <- rownames(pval)[rowSums(pval <= cutoff) != 0]
  pval <- pval[sig_inter, ]

  if (cluster_proteins) {
    ord_proteins <- rownames(pval)[hclust(dist(pval))$order]
  }

  if (cluster_cells) {
    ord_cells <- colnames(pval)[hclust(dist(t(pval)))$order]
  }

  d.pval <- d.pval |> filter(.data[["interacting_pair"]] %in% sig_inter)
  d.mean <- d.mean |> filter(.data[["interacting_pair"]] %in% sig_inter)

  d.pval <- d.pval |> pivot_longer(cells, names_to="cells", values_to="p.value")
  d.mean <- d.mean |> pivot_longer(cells, names_to="cells", values_to="mean")

  d <- left_join(d.pval, d.mean, by = c("interacting_pair", "cells"))

  if (cluster_proteins)
    d$interacting_pair <- factor(d$interacting_pair, ord_proteins)

  if (cluster_cells)
    d$cells <- factor(d$cells, ord_cells)

  ggplot(d, aes(interacting_pair, cells, size=-log10(p.value + 1e-3), color=log2(mean))) +
    geom_point() +
    scale_color_distiller(palette="RdYlBu") +
    guides(size=guide_legend("-log10(p.value)")) +
    Seurat::RotatedAxis()
}

#' plot_cpdb_heatmap
#'
#' Plot cellphonedb results as a heatmap of cell types.
#'
#' @param path path to the folder containing cellphonedb results.
#' @param cutoff cutoff for pvalues.
#'
#' @details The score is calculated as the number of significant interactions (based on cutoff) per cell pair.
#'
#' @export
#'
plot_cpdb_heatmap <- function(x, cutoff=1e-3) {
  d <- calculate_cpdb_cell_graph(x, cutoff=cutoff)
  ggplot(d, aes(.data[["cell_b"]], .data[["cell_a"]], fill=.data[["score"]])) +
    geom_tile() +
    scale_fill_distiller(palette="Spectral") +
    Seurat::RotatedAxis()
}

calculate_cpdb_cell_graph <- function(x, cutoff=1e-3) {
  res <- colSums(x[["m_pval"]] <= cutoff)
  d <- data.frame(cells=names(res), score=res)
  d <- d |> separate(col="cells", into=c("cell_a", "cell_b"), sep="\\|")
  m <- d |> pivot_wider(names_from="cell_b", values_from="score", values_fill=0) |> column_to_rownames("cell_a")
  o <- hclust(dist(m))$order
  cell_order <- rownames(m)[o]
  d |> mutate(cell_a=factor(.data[["cell_a"]], levels=cell_order)) |>
    mutate(cell_b=factor(.data[["cell_b"]], levels=cell_order))
}

#' @export
calculate_cpdb_cellxcell_matrix <- function(x) {
  calculate_cpdb_cell_graph(x) |> pivot_wider(names_from="cell_b", values_from="score") |> column_to_rownames("cell_a") |> as.matrix()
}

#' @export
plot_cpdb_cellxcell <- function(x, name="scores", ...) {
    m <- calculate_cpdb_cellxcell_matrix(x)
    Heatmap(m, name=name, ...)
}

#' @export
plot_cpdb_cellxpair <- function(x, cells=NULL, pval.cutoff=0.01, mean.cutoff=1, scale=FALSE, ...) {
  m_mean <- as.matrix(x$m_mean)
  rownames(m_mean) <- x$proteins

  m_pval <- as.matrix(x$m_pval)

  if (is.null(cells))
    cells <- x$cells

  if (scale)
    s_mean <- t(scale(t(m_mean)))
  else
    s_mean <- m_mean

  sel.pval <- rowSums(m_pval[, cells, drop=FALSE] < pval.cutoff) > 0
  sel.mean <- rowSums(m_mean[, cells, drop=FALSE] > mean.cutoff) > 0

  Heatmap(s_mean[sel.pval & sel.mean, cells], show_column_names=TRUE, show_row_names=TRUE, name="score", ...)
}
