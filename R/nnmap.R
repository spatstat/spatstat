#
#  nnmap.R
#
#    k-th nearest neighbour of each pixel
#
#  $Revision: 1.2 $  $Date: 2013/05/26 07:39:46 $
#

nnmap <- function(X, k=1, what = c("dist", "which"), ...,
                  W=as.owin(X),
                  is.sorted.X=FALSE,
                  sortby=c("range", "var", "x", "y")) {
  stopifnot(is.ppp(X))
  sortby <- match.arg(sortby)
  outputarray <- resolve.1.default("outputarray", ..., outputarray=FALSE)
  
  W <- as.owin(W)
  huge <- 1.1 * diameter(bounding.box(as.rectangle(X), as.rectangle(W)))
  
  what   <- match.arg(what, choices=c("dist", "which"), several.ok=TRUE)
  want.dist  <- "dist" %in% what 
  want.which <- "which" %in% what
  want.both  <- want.dist && want.which

  if(!missing(k)) {
    # k can be a single integer or an integer vector
    if(length(k) == 0)
      stop("k is an empty vector")
    else if(length(k) == 1) {
      if(k != round(k) || k <= 0)
        stop("k is not a positive integer")
    } else {
      if(any(k != round(k)) || any(k <= 0))
        stop(paste("some entries of the vector",
                   sQuote("k"), "are not positive integers"))
    }
  }
  k <- as.integer(k)
  kmax <- max(k)
  nk <- length(k)

  # set up pixel array
  M <- do.call.matched("as.mask", append(list(w=W), list(...)))
  Mdim <- M$dim
  nxcol <- Mdim[2]
  nyrow <- Mdim[1]
  npixel <- nxcol * nyrow
  
  nX <- npoints(X)
  if(nX == 0) {
    # trivial - avoid potential problems in C code
    NND <- if(want.dist) array(Inf, dim=c(nk, Mdim)) else 0
    NNW <- if(want.which) array(NA_integer_, dim=c(nk, Mdim)) else 0
  } else {
    # usual case 
    if(is.sorted.X && !(sortby %in% c("x", "y")))
      stop(paste("If data are already sorted,",
                 "the sorting coordinate must be specified explicitly",
                 "using sortby = \"x\" or \"y\""))

    # decide whether to sort on x or y coordinate
    switch(sortby,
           range = {
             s <- sidelengths(as.rectangle(X))
             sortby.y <- (s[1] < s[2])
           },
           var = {
             sortby.y <- (var(X$x) < var(X$y))
           },
           x={ sortby.y <- FALSE},
           y={ sortby.y <- TRUE}
           )

    # The C code expects points to be sorted by x coordinate.
    if(sortby.y) {
      oldM <- M
      X <- flipxy(X)
      W <- flipxy(W)
      M <- flipxy(M)
      Mdim <- M$dim
    }
    xx <- X$x
    yy <- X$y
    # sort only if needed
    if(!is.sorted.X){
      oX <- fave.order(xx)
      xx <- xx[oX]
      yy <- yy[oX]
    }

    # number of neighbours that are well-defined
    kmaxcalc <- min(nX, kmax)

    # call C code
    Cfun <- paste0("knnG",
                   if(want.dist) "d" else "",
                   if(want.which) "w" else "")
    nndv <- if(want.dist) numeric(npixel * kmaxcalc) else numeric(1)
    nnwh <- if(want.which) integer(npixel * kmaxcalc) else integer(1)

    DUP <- spatstat.options("dupC")
   
    zz <- .C(Cfun,
             nx = as.integer(nxcol),
             x0 = as.double(M$xcol[1]),
             xstep = as.double(M$xstep),
             ny = as.integer(nyrow),
             y0 = as.double(M$yrow[1]),
             ystep = as.double(M$ystep),
             np = as.integer(nX),
             xp = as.double(xx),
             yp = as.double(yy),
             kmax = as.integer(kmaxcalc),
             nnd = as.double(nndv),
             nnwhich = as.integer(nnwh),
             huge = as.double(huge),
             DUP = DUP,
             PACKAGE = "spatstat")

    # extract results
    nnW <- zz$nnwhich
    nnD <- zz$nnd
    # map index 0 to NA
    if(want.which && any(uhoh <- (nnW == 0))) {
      nnW[uhoh] <- NA
      if(want.dist) nnD[uhoh] <- Inf
    }
    # reinterpret indices in original ordering
    if(!is.sorted.X) nnW <- oX[nnW]
  
    # reform as arrays 
    NND <- if(want.dist) array(nnD, dim=c(kmaxcalc, Mdim)) else 0
    NNW <- if(want.which) array(nnW, dim=c(kmaxcalc, Mdim)) else 0
    if(sortby.y) {
      # flip x and y back again
      if(want.dist) NND <- aperm(NND, c(1, 3, 2))
      if(want.which) NNW <- aperm(NNW, c(1, 3, 2))
      M <- oldM
      Mdim <- dim(M)
    }
    
    # the return value should correspond to the original vector k
    if(kmax > kmaxcalc) {
      # pad with NA / Inf
      if(want.dist) {
        NNDcalc <- NND
        NND <- array(Inf, dim=c(kmax, Mdim))
        NND[1:kmaxcalc, , ] <- NNDcalc
      }
      if(want.which) {
        NNWcalc <- NNW
        NNW <- array(NA_integer_, dim=c(kmax, Mdim))
        NNW[1:kmaxcalc, , ] <- NNWcalc
      }
    }
    if(length(k) < kmax) {
      # select only the specified planes
      if(want.dist)
        NND <- NND[k, , , drop=FALSE]
      if(want.which)
        NNW <- NNW[k, , , drop=FALSE]
    }
  }

  # secret backdoor
  if(outputarray) {
    # return result as an array or pair of arrays
    result <- if(want.both) { list(dist=NND, which=NNW) } else
              if(want.dist) NND else NNW
    attr(result, "pixarea") <- with(M, xstep * ystep)
    return(result)
  }

  # format result as a list of images
  result <- list()
  if(want.dist) {
    dlist <- list()
    for(i in 1:nk)
      dlist[[i]] <- as.im(NND[i,,], M)
    names(dlist) <- k
    result[["dist"]] <- if(nk > 1) dlist else dlist[[1]]
  }
  if(want.which) {
    wlist <- list()
    for(i in 1:nk)
      wlist[[i]] <- as.im(NNW[i,,], M)
    names(wlist) <- k
    result[["which"]] <- if(nk > 1) wlist else wlist[[1]]
  }
  if(!want.both) result <- result[[1]]
  return(result)
}
