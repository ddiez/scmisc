#' compute_percentage_variance
#'
#' @param x and object with PCA coordintates
#' @param reduction name of the reduction to use with Seurat objects.
#'
#' @export
#'
compute_percentage_variance <- function(x, reduction = "pca") {
  UseMethod("compute_percentage_variance")
}


#' @rdname compute_percentage_variance
#' @export
compute_percentage_variance.Seurat <- function(x, reduction = "pca") {
  emb <- Embeddings(x, reduction = reduction)
  stdev <- Stdev(x, reduction = reduction)


  d <- tibble(dims = seq_len(ncol(emb)), stdev = stdev, variance = stdev ^ 2) |>
    add_tally(.data$variance, name = "total") |>
    mutate(percentage = 100 * .data$variance / .data$total)
}

