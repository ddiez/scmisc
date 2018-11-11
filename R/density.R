#' compute_density
#'
#' @param x x coordinate of data
#' @param y y coordinate of data
#' @param n Number of grid points in each direction. Can be scalar or a length-2 integer vector.
#' @param coord.name name of coordinates for some methods.
#' @param ... additional arguments passed down to methods.
#'
#' @note I took this implementation from http://slowkow.com/notes/ggplot2-color-by-density/
#' @export
compute_density <- function(x, ...) {
  UseMethod("compute_density")
}

#' @rdname compute_density
#' @export
compute_density.default <- function(x, y, n = 100, ...) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  dens$z[ii]
}

#' @rdname compute_density
#' @export
compute_density.SingleCellExperiment <- function(x, coord.name = "TSNE", ...) {
  y <- get_coord(x, coord.name = coord.name, add.cols = FALSE)
  compute_density(y[, 1], y[, 2])
}

#' @rdname compute_density
#' @export
compute_density.CellDataSet <- function(x, coord.name = "A", ...) {
  y <- get_coord(x, coord.name = coord.name, add.cols = FALSE)
  compute_density(y[, 1], y[, 2])
}

#' @rdname compute_density
#' @export
cell_density <- function(x, ...) {
  compute_density(x, ...)
}
