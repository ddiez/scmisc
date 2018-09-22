#' get_density
#'
#' @param x x coordinate of data
#' @param y y coordinate of data
#' @param n Number of grid points in each direction. Can be scalar or a length-2 integer vector.
#'
#' @note I took this implementation from http://slowkow.com/notes/ggplot2-color-by-density/
#' @export
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  dens$z[ii]
}
