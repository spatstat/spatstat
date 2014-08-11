#
#  cut.ppp.R
#
#  cut method for ppp objects
#
#  $Revision: 1.12 $   $Date: 2013/04/25 06:37:43 $
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
    } else if((length(z) == 1) && (z %in% colnames(marks(x)))) {
      # interpret as the name of a column of marks
      zname <- z
      m <- marks(x)
      z <- m[, zname]
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

  if(is.tess(z)) {
    switch(z$type,
           rect={
             jx <- findInterval(x$x, z$xgrid, rightmost.closed=TRUE)
             iy <- findInterval(x$y, z$ygrid, rightmost.closed=TRUE)
             nrows    <- length(z$ygrid) - 1
             ncols <- length(z$xgrid) - 1
             jcol <- jx
             irow <- nrows - iy + 1
             ktile <- jcol + ncols * (irow - 1)
             m <- factor(ktile, levels=seq_len(nrows*ncols))
             ij <- expand.grid(j=seq_len(ncols),i=seq_len(nrows))
             levels(m) <- paste("Tile row ", ij$i, ", col ", ij$j, sep="")
           },
           tiled={
             todo <- seq_len(npoints(x))
             nt <- length(z$tiles)
             m <- integer(x$n)
             for(i in 1:nt) {
               ti <- z$tiles[[i]]
               hit <- inside.owin(x$x[todo], x$y[todo], ti)
               if(any(hit)) {
                 m[todo[hit]] <- i
                 todo <- todo[!hit]
               }
               if(length(todo) == 0)
                 break
             }
             m[m == 0] <- NA
             nama <- names(z$tiles)
             lev <- seq_len(nt)
             lab <- if(!is.null(nama) && all(nzchar(nama))) nama else paste("Tile", lev)
             m <- factor(m, levels=lev, labels=lab)
           },
           image={
             zim <- z$image
             m <- factor(zim[x, drop=FALSE], levels=levels(zim))
           }
           )
    marks(x) <- m
    return(x)
  }
  stop("Format of z not understood")
} 

