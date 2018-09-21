#' fortify.DataFrame
#'
#' fortify method for DataFrame objects.
#'
#' @param model model or other R object to convert to data frame
#' @param data original dataset, if needed
#' @param ... other arguments passed to methods
#'
#' @export
#'
fortify.DataFrame <- function(model, data, ...) {
  as.data.frame(model)
}
