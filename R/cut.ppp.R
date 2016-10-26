#
#  cut.ppp.R
#
#  cut method for ppp objects
#
#  $Revision: 1.15 $   $Date: 2016/10/26 09:29:57 $
#

cut.ppp <- function(x, z=marks(x), ...) {
  x <- as.ppp(x)
  if(missing(z) || is.null(z)) {
    z <- marks(x, dfok=TRUE)
    if(is.null(z))
      stop("x has no marks to cut")
  }
  if(is.character(z)) {
    if(length(z) == npoints(x)) {
      # interpret as a factor
      z <- factor(z)
    } else if((length(z) == 1) && (z %in% colnames(df <- as.data.frame(x)))) {
      # interpret as the name of a column of marks or coordinates
      z <- df[, z]
    } else stop("format of argument z not understood") 
  }
  if(is.factor(z) || is.vector(z)) {
    stopifnot(length(z) == npoints(x))
    g <- if(is.factor(z)) z else if(is.numeric(z)) cut(z, ...) else factor(z)
    marks(x) <- g
    return(x)
  }
  if(is.data.frame(z) || is.matrix(z)) {
    stopifnot(nrow(z) == npoints(x))
    # take first column 
    z <- z[,1]
    g <- if(is.numeric(z)) cut(z, ...) else factor(z)
    marks(x) <- g
    return(x)
  }
  if(is.im(z)) 
    return(cut(x, z[x, drop=FALSE], ...))

  if(is.owin(z)) {
    marks(x) <- factor(inside.owin(x$x, x$y, z), levels=c(FALSE, TRUE))
    return(x)
  }
  
  if(is.tess(z)) {
    marks(x) <- tileindex(x$x, x$y, z)
    return(x)
  }

  stop("Format of z not understood")
} 

