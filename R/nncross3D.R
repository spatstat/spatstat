#
#   nncross3D.R
#
#    $Revision: 1.2 $  $Date: 2013/06/29 03:30:45 $
#
#  Copyright (C) Adrian Baddeley, Jens Oehlschlaegel and Rolf Turner 2000-2013
#  Licence: GNU Public Licence >= 2

nncross.pp3 <- function(X, Y, iX=NULL, iY=NULL,
                    what = c("dist", "which"),
                    ...,
                    k = 1,
                    sortby=c("range", "var", "x", "y", "z"),
                    is.sorted.X = FALSE,
                    is.sorted.Y = FALSE) {
  stopifnot(is.pp3(Y))
  sortby <- match.arg(sortby)
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

  # trivial cases
  nX <- npoints(X)
  nY <- nobjects(Y)
  # deal with null cases
  if(nX == 0)
    return(as.data.frame(list(dist=matrix(0, nrow=0, ncol=nk),
                which=matrix(0L, nrow=0, ncol=nk))[what]))
  if(nY == 0)
    return(as.data.frame(list(dist=matrix(Inf, nrow=nX, ncol=nk),
                             which=matrix(NA, nrow=nX, ncol=nk))[what]))
  
  if(is.null(iX) != is.null(iY))
    stop("If one of iX, iY is given, then both must be given")
  exclude <- (!is.null(iX) || !is.null(iY))
  if(exclude) {
    stopifnot(is.integer(iX) && is.integer(iY))
    if(length(iX) != nX)
      stop("length of iX does not match the number of points in X")
    if(length(iY) != nY)
      stop("length of iY does not match the number of points in Y")
  }

  if((is.sorted.X || is.sorted.Y) && !(sortby %in% c("x", "y", "z")))
     stop(paste("If data are already sorted,",
                "the sorting coordinate must be specified explicitly",
                "using sortby = \"x\" or \"y\" or \"z\""))

  # decide which coordinate to sort on
  switch(sortby,
         range = {
           s <- sidelengths(as.box3(Y))
           sortcoord <- c("x", "y", "z")[which.min(s)]
         },
         var = {
           v <- apply(coords(Y), 2, var)
           sortcoord <- c("x", "y", "z")[which.min(v)]           
         },
         x={ sortcoord <- "x" },
         y={ sortcoord <- "y" },
         z={ sortcoord <- "z" }
         )

  # The C code expects points to be sorted by z coordinate.
  XX <- coords(X)
  YY <- coords(Y)
  switch(sortcoord,
         x = {
           # rotate x axis to z axis
           XX <- XX[, c(3,2,1)]
           YY <- YY[, c(3,2,1)]
         },
         y = {
           # rotate y axis to z axis
           XX <- XX[, c(3,1,2)]
           YY <- YY[, c(3,1,2)]
         },
         z = { })

  # sort only if needed
  if(!is.sorted.X){
    oX <- fave.order(XX[,3])
    XX[oX, , drop=FALSE]
    if(exclude) iX <- iX[oX]
  }
  if (!is.sorted.Y){
    oY <- fave.order(YY[,3])
    YY <- YY[oY, , drop=FALSE]
    if(exclude) iY <- iY[oY]
  }

  # number of neighbours that are well-defined
  kmaxcalc <- min(nY, kmax)
  
  if(kmaxcalc == 1) {
    # ............... single nearest neighbour ..................
    # call C code
    Cfun <- paste0("nnX",
                   if(exclude) "E" else "",
                   if(want.dist) "d" else "",
                   if(want.which) "w" else "",
                   "3D")
    nndv <- if(want.dist) numeric(nX) else numeric(1)
    nnwh <- if(want.which) integer(nX) else integer(1)
    if(!exclude) iX <- iY <- integer(1)
    DUP <- spatstat.options("dupC")
    huge <- 1.1 * diameter(bounding.box3(as.box3(X),as.box3(Y)))
  
    z <- .C(Cfun,
            n1=as.integer(nX),
            x1=as.double(XX[,1]),
            y1=as.double(XX[,2]),
            z1=as.double(XX[,3]),
            id1=as.integer(iX),
            n2=as.integer(nY),
            x2=as.double(YY[,1]),
            y2=as.double(YY[,2]),
            z2=as.double(YY[,3]),
            id2=as.integer(iY),
            nnd=as.double(nndv),
            nnwhich=as.integer(nnwh),
            huge=as.double(huge),
            DUP=DUP,
            PACKAGE="spatstat")

    if(want.which) {
      nnwcode <- z$nnwhich #sic. C code now increments by 1
      if(any(uhoh <- (nnwcode == 0))) {
        warning("NA's produced in nncross()$which")
        nnwcode[uhoh] <- NA
      }
    }
  
    # reinterpret in original ordering
    if(is.sorted.X){
      if(want.dist) nndv <- z$nnd
      if(want.which) nnwh <- if(is.sorted.Y) nnwcode else oY[nnwcode]
    } else {
      if(want.dist) nndv[oX] <- z$nnd
      if(want.which) nnwh[oX] <- if(is.sorted.Y) nnwcode else oY[nnwcode]
    }

    if(want.both) return(data.frame(dist=nndv, which=nnwh))
    return(if(want.dist) nndv else nnwh)

  } else {
    # ............... k nearest neighbours ..................
    # call C code
    Cfun <- paste0("knnX",
                  if(exclude) "E" else "",
                   if(want.dist) "d" else "",
                   if(want.which) "w" else "",
                   "3D")
    nndv <- if(want.dist) numeric(nX * kmaxcalc) else numeric(1)
    nnwh <- if(want.which) integer(nX * kmaxcalc) else integer(1)
    if(!exclude) iX <- iY <- integer(1)
    DUP <- spatstat.options("dupC")
    huge <- 1.1 * diameter(bounding.box3(as.box3(X),as.box3(Y)))
  
    z <- .C(Cfun,
            n1=as.integer(nX),
            x1=as.double(XX[,1]),
            y1=as.double(XX[,2]),
            z1=as.double(XX[,3]),
            id1=as.integer(iX),
            n2=as.integer(nY),
            x2=as.double(YY[,1]),
            y2=as.double(YY[,2]),
            z2=as.double(YY[,3]),
            id2=as.integer(iY),
            kmax=as.integer(kmaxcalc),
            nnd=as.double(nndv),
            nnwhich=as.integer(nnwh),
            huge=as.double(huge),
            DUP=DUP,
            PACKAGE="spatstat")

    # extract results
    nnD <- z$nnd
    nnW <- z$nnwhich
    # map 0 to NA
    if(want.which && any(uhoh <- (nnW == 0))) {
      nnW[uhoh] <- NA
      if(want.dist) nnD[uhoh] <- Inf
    }
    # reinterpret indices in original ordering
    if(!is.sorted.Y) nnW <- oY[nnW]
    # reform as matrices
    NND <- if(want.dist) matrix(nnD, nrow=nX, ncol=kmaxcalc, byrow=TRUE) else 0
    NNW <- if(want.which) matrix(nnW, nrow=nX, ncol=kmaxcalc, byrow=TRUE) else 0
    if(!is.sorted.X){
      # rearrange rows to correspond to original ordering of points
      if(want.dist) NND[oX, ] <- NND
      if(want.which) NNW[oX, ] <- NNW
    }
    # the return value should correspond to the original vector k
    if(kmax > kmaxcalc) {
      # add columns of NA / Inf
      kextra <- kmax - kmaxcalc
      if(want.dist)
        NND <- cbind(NND, matrix(Inf, nrow=nX, ncol=kextra))
      if(want.which)
        NNW <- cbind(NNW, matrix(NA_integer_, nrow=nX, ncol=kextra))
    }
    if(length(k) < kmax) {
      # select only the specified columns
      if(want.dist)
        NND <- NND[, k, drop=TRUE]
      if(want.which)
        NNW <- NNW[, k, drop=TRUE]
    }

    result <- as.data.frame(list(dist=NND, which=NNW)[what])
    if(ncol(result) == 1)
      result <- result[, , drop=TRUE]
    return(result)
  }
}

