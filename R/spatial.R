#' plot_spatial
#' @export
plot_spatial <- function(x, features, size=.5, rotate=FALSE) {
  UseMethod("plot_spatial")
}

#' @rdname plot_spatial
#' @export
plot_spatial.Seurat <- function(x, features, size=.5, rotate=FALSE) {
  images <- SeuratObject::Images(x)
  p <- vector("list", length(images) * length(features))
  k <- 1
  exprs <- SeuratObject::GetAssayData(x, slot="data")[features, , drop=FALSE]
  #maxlim <- apply(exprs, 1, max)
  #names(maxlim) <- features
  for (image in images) {
    coord <- SeuratObject::GetTissueCoordinates(x, image=image)

    sel <- grep(paste0(image, "_"), colnames(exprs))
    d <- cbind(coord, Matrix::t(exprs[, sel, drop=FALSE]))

    if (rotate) {
      d$x <- -d$imagerow
      d$y <- -d$imagecol
    }
    else {
      d$x <- d$imagecol
      d$y <- -d$imagerow
    }

    for (feature in features) {
      p[[k]] <- ggplot(d, aes(.data[["x"]], .data[["y"]], fill=.data[[feature]])) +
        geom_point(shape=21, size=size, stroke=.1) +
        scale_fill_gradientn(colors=Seurat:::SpatialColors(1000)) + #, limits=c(NA, maxlim[feature])) +
        labs(title=image, subtitle=feature) +
        theme(axis.line=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank())
      k <- k + 1
    }
  }
  wrap_plots(p, byrow=FALSE, ncol=length(images)) & theme(aspect.ratio=1)
}
