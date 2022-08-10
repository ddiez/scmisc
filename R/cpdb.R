#' plot_cpdb_dotplot
#'
#' Plot cellphonedb results as a dotplot.
#'
#' @param path path to the folder containing cellphonedb results.
#' @param cutoff cutoff for pvalues.
#' @param filter filter for cell pairs. Can be a regular expression.
#'
#' @export
#'
plot_cpdb_dotplot <- function(path, cutoff=1e-3, filter=NULL) {
  cols <- c("id_cp_interaction", "interacting_pair", "partner_a", "partner_b", "gene_a", "gene_b", "secreted", "receptor_a", "receptor_b", "annotation_strategy", "is_integrin")
  d.pval <- read_tsv(file.path(path, "pvalues.txt"), show_col_types=FALSE)
  d.mean <- read_tsv(file.path(path, "means.txt"), show_col_types=FALSE)

  cells <- setdiff(colnames(d.pval), cols)

  if (!is.null(filter))
    cells <- grep(filter, cells, value=TRUE)

  d.pval <- d.pval[, c("interacting_pair", cells)]
  d.mean <- d.mean[, c("interacting_pair", cells)]

  d.pval[["interacting_pair"]] <- make.names(d.pval[["interacting_pair"]], unique=TRUE)
  d.mean[["interacting_pair"]] <- make.names(d.mean[["interacting_pair"]], unique=TRUE)

  pval <- d.pval |> column_to_rownames("interacting_pair") |> as.matrix()
  sig_inter <- rownames(pval)[rowSums(pval < cutoff) != 0]

  d.pval <- d.pval |> filter(.data[["interacting_pair"]] %in% sig_inter)
  d.mean <- d.mean |> filter(.data[["interacting_pair"]] %in% sig_inter)

  d.pval <- d.pval |> pivot_longer(cells, names_to="cells", values_to="p.value")
  d.mean <- d.mean |> pivot_longer(cells, names_to="cells", values_to="mean")

  d <- left_join(d.pval, d.mean, by = c("interacting_pair", "cells"))
  ggplot(d, aes(interacting_pair, cells, size=-log10(p.value + 1e-3), color=log2(mean))) +
    geom_point() +
    scale_color_distiller(palette="RdYlBu") +
    guides(size=guide_legend("-log10(p.value)")) +
    Seurat::RotatedAxis()
}
