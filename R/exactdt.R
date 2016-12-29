#
#	exactdt.S
#	S function exactdt() for exact distance transform
#
#	$Revision: 4.16 $	$Date: 2014/10/24 00:22:30 $
#

exactdt <- local({

  die <- function(why) { stop(paste("ppp object format corrupted:", why)) }

  exactdt <- function(X, ...) {
    verifyclass(X, "ppp")
    w <- X$window
    if(spatstat.options("exactdt.checks.data")) {
      ## check validity of ppp structure 
      bb <- as.rectangle(w)
      xr <- bb$xrange
      yr <- bb$yrange
      rx <- range(X$x)
      ry <- range(X$y)
      if(rx[1L] < xr[1L] || rx[2L] > xr[2L]) die("x-coordinates out of bounds")
      if(ry[1L] < yr[1L] || ry[2L] > yr[2L]) die("y-coordinates out of bounds")
      if(length(X$x) != length(X$y))
        die("x and y vectors have different length")
      if(length(X$x) != X$n) die("length of x,y vectors does not match n")
    }
    w <- as.mask(w, ...)
    ## dimensions of result
    nr <- w$dim[1L]
    nc <- w$dim[2L]
    ## margins in C array 
    mr <- 2
    mc <- 2
    ## full dimensions of allocated storage
    Nnr <- nr + 2 * mr
    Nnc <- nc + 2 * mc
    N <- Nnr * Nnc
    ## output rows & columns (R indexing)
    rmin <- mr + 1
    rmax <- Nnr - mr
    cmin <- mc + 1
    cmax <- Nnc - mc
    ## go
    res <- .C("exact_dt_R",
              as.double(X$x),
              as.double(X$y),
              as.integer(X$n),
              as.double(w$xrange[1L]),
              as.double(w$yrange[1L]),
              as.double(w$xrange[2L]),
              as.double(w$yrange[2L]),
              nr = as.integer(nr),
              nc = as.integer(nc),
              mr = as.integer(mr),
              mc = as.integer(mc),
              distances = as.double(double(N)),
              indices = as.integer(integer(N)),
              boundary = as.double(double(N)))
    ## extract 
    dist <- matrix(res$distances,
                   ncol=Nnc, nrow=Nnr, byrow = TRUE)[rmin:rmax, cmin:cmax]
    inde <- matrix(res$indices,
                   ncol=Nnc, nrow=Nnr, byrow = TRUE)[rmin:rmax, cmin:cmax]
    bdry <- matrix(res$boundary,
                   ncol=Nnc, nrow=Nnr, byrow = TRUE)[rmin:rmax, cmin:cmax]
    ## convert index from C to R indexing
    inde <- inde + 1L
    return(list(d = dist, i = inde, b = bdry, w=w))
  }

  exactdt
})

