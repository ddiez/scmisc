#' @rdname plot_spatial
#' @export
plot_spatial <- function(x, features, size=.5) {
  UseMethod("plot_spatial")
}

#' @rdname plot_spatial
#' @export
plot_spatial.Seurat <- function(x, features, size=.5) {
  images <- Seurat::Images(x)
  p <- vector("list", length(images) * length(features))
  k <- 1
  exprs <- GetAssayData(x, slot="data")[features, ]
  #maxlim <- apply(exprs, 1, max)
  #names(maxlim) <- features
  for (image in images) {
    coord <- Seurat::GetTissueCoordinates(x, image=image)

    sel <- grep(paste0(image, "_"), colnames(exprs))
    d <- cbind(coord, Matrix::t(exprs[, sel]))

    for (feature in features) {
      p[[k]] <- ggplot(d, aes(-imagerow, imagecol, fill=.data[[feature]])) +
        geom_point(shape=21, size=size, stroke=.1) +
        coord_flip() +
        scale_fill_gradientn(colors=Seurat:::SpatialColors(1000)) + #, limits=c(NA, maxlim[feature])) +
        labs(title=image, subtitle=feature) +
        NoAxes()
      k <- k + 1
    }
  }
  wrap_plots(p, byrow=FALSE, ncol=length(images)) + theme(aspect.ratio=1)
}
