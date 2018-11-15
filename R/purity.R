#' compute_purity
#'
#' @param x object to compute purity on.
#' @param col.x name of column.
#' @param col.y name of column.
#' @param ...
#'
#' @export
compute_purity <- function(x, ...) {
  UseMethod("compute_purity")
}

#' @rdname compute_purity
#' @export
compute_purity.DataFrame <- function(x, ...) {
  compute_purity(as.data.frame(x), ...)
}

#' @rdname compute_purity
#' @export
compute_purity.data.frame <- function(x, col.x, col.y, ...) {
  if (missing(col.x)) stop("col.x is missing.")
  if (missing(col.y)) stop("col.y is missing.")

  x %>% select_(col.x, col.y) %>%
    group_by_(col.y) %>%
    mutate(total = n()) %>%
    group_by_(col.x, col.y) %>%
    summarize(count = n(), purity = count / first(.data$total))
}

