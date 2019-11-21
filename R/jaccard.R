#' compute_jaccard
#'
#' @param x object to compute jaccard index on.
#' @param col.x name of column.
#' @param col.y name of column.
#' @param ... arguments passed down to methods.
#'
#' @export
compute_jaccard <- function(x, ...) {
  UseMethod("compute_jaccard")
}

#' @rdname compute_jaccard
#' @export
compute_jaccard.SingleCellExperiment <- function(x, ...) {
  compute_jaccard(colData(x), ...)
}

#' @rdname compute_jaccard
#' @export
compute_jaccard.DataFrame <- function(x, ...) {
  compute_jaccard(as.data.frame(x), ...)
}

#' @rdname compute_jaccard
#' @export
compute_jaccard.data.frame <- function(x, col.x, col.y, ...) {
  if (missing(col.x)) stop("col.x is missing.")
  if (missing(col.y)) stop("col.y is missing.")

  d <- x %>% select_at(.vars = c(col.x, col.y))
  d <- d %>% group_by_at(.vars = c(col.x, col.y)) %>% dplyr::count() %>% ungroup()
  d <- d %>% group_by_at(.vars = col.y) %>%
    mutate(total_y = sum(.data$n)) %>%
    group_by_at(.vars = col.x) %>%
    mutate(total_x = sum(.data$n))
  d %>%
    mutate(jaccard = 2 * .data$n / (.data$total_x + .data$total_y)) %>%
    select_at(.vars = c(col.x, col.y, count =  "n", "jaccard"))
}
