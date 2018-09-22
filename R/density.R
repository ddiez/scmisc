#' get_density
#'
#' @param x x coordinate of data
#' @param y y coordinate of data
#' @param n Number of grid points in each direction. Can be scalar or a length-2 integer vector.
#'
#' @note I took this implementation from http://slowkow.com/notes/ggplot2-color-by-density/
#' @export
get_density <- function(x, ...) {
  UseMethod("get_density")
}

#' @rdname get_density
#' @export
get_density.default <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  dens$z[ii]
}

#' @rdname get_density
#' @export
get_density.SingleCellExperiment <- function(x, name = "tSNE") {
  y <- get_coord(x, name = name, annotate = FALSE)
  get_density(y[, 1], y[, 2])
}

#' @rdname get_density
#' @export
get_density.CellDataSet <- function(x, name = "A") {
  y <- get_coord(x, name = name, annotate = FALSE)
  get_density(y[, 1], y[, 2])
}
