#
#    pspcross.R
#
#    Intersections of line segments
#    
#    $Revision: 1.20 $   $Date: 2015/10/21 09:06:57 $
#
#
crossing.psp <- function(A,B,fatal=TRUE) {
  verifyclass(A, "psp")
  verifyclass(B, "psp")
  
  # first check for intersection of windows
  ABW <- intersect.owin(A$window, B$window, fatal=fatal)
  if(is.null(ABW)) 
    return(NULL)
  
  eps <- .Machine$double.eps

  na <- A$n
  eA <- A$ends
  x0a <- eA$x0
  y0a <- eA$y0
  dxa <- eA$x1 - eA$x0
  dya <- eA$y1 - eA$y0

  nb <- B$n
  eB <- B$ends
  x0b <- eB$x0
  y0b <- eB$y0
  dxb <- eB$x1 - eB$x0
  dyb <- eB$y1 - eB$y0

  useCall <- spatstat.options("crossing.psp.useCall")
  if(!useCall) {
    # old C routine
    out <- .C("xysegint",
              na=as.integer(na),
              x0a=as.double(x0a),
              y0a=as.double(y0a),
              dxa=as.double(dxa),
              dya=as.double(dya), 
              nb=as.integer(nb),
              x0b=as.double(x0b),
              y0b=as.double(y0b),
              dxb=as.double(dxb),
              dyb=as.double(dyb), 
              eps=as.double(eps),
              xx=as.double(numeric(na * nb)),
              yy=as.double(numeric(na * nb)),
              ta=as.double(numeric(na * nb)),
              tb=as.double(numeric(na * nb)),
              ok=as.integer(integer(na * nb)))
    
    ok <- (matrix(out$ok, na, nb) != 0)
    xx <- matrix(out$xx, na, nb)
    yy <- matrix(out$yy, na, nb)
    xx <- as.vector(xx[ok])
    yy <- as.vector(yy[ok])
  } else {
    # new
    storage.mode(x0a) <- storage.mode(y0a) <- "double"
    storage.mode(dxa) <- storage.mode(dya) <- "double"
    storage.mode(x0b) <- storage.mode(y0b) <- "double"
    storage.mode(dxb) <- storage.mode(dyb) <- "double"
    storage.mode(eps) <- "double"
    out <- .Call("Cxysegint",
                 x0a, 
                 y0a, 
                 dxa, 
                 dya, 
                 x0b, 
                 y0b, 
                 dxb, 
                 dyb, 
    	         eps)
#                 PACKAGE="spatstat")
    xx <- out[[5]]
    yy <- out[[6]]
  }
  result <- ppp(xx, yy, window=ABW, check=FALSE)
  return(result)
}

test.crossing.psp <- function(A,B) {
  # return logical matrix specifying whether A[i] and B[j] cross
  verifyclass(A, "psp")
  verifyclass(B, "psp")
  eps <- .Machine$double.eps

  na <- A$n
  eA <- A$ends
  x0a <- eA$x0
  y0a <- eA$y0
  dxa <- eA$x1 - eA$x0
  dya <- eA$y1 - eA$y0

  nb <- B$n
  eB <- B$ends
  x0b <- eB$x0
  y0b <- eB$y0
  dxb <- eB$x1 - eB$x0
  dyb <- eB$y1 - eB$y0

  out <- .C("xysi",
            na=as.integer(na),
            x0a=as.double(x0a),
            y0a=as.double(y0a),
            dxa=as.double(dxa),
            dya=as.double(dya), 
            nb=as.integer(nb),
            x0b=as.double(x0b),
            y0b=as.double(y0b),
            dxb=as.double(dxb),
            dyb=as.double(dyb), 
            eps=as.double(eps),
            ok=as.integer(integer(na * nb)))

  hit <- (matrix(out$ok, na, nb) != 0)
  return(hit)
}

anycrossing.psp <- function(A,B) {
  # equivalent to: any(test.crossing.psp(A,B))
  # Test whether two psp objects have at least one crossing point
  verifyclass(A, "psp")
  verifyclass(B, "psp")
  eps <- .Machine$double.eps

  na <- A$n
  eA <- A$ends
  x0a <- eA$x0
  y0a <- eA$y0
  dxa <- eA$x1 - eA$x0
  dya <- eA$y1 - eA$y0

  nb <- B$n
  eB <- B$ends
  x0b <- eB$x0
  y0b <- eB$y0
  dxb <- eB$x1 - eB$x0
  dyb <- eB$y1 - eB$y0

  out <- .C("xysiANY",
            na=as.integer(na),
            x0a=as.double(x0a),
            y0a=as.double(y0a),
            dxa=as.double(dxa),
            dya=as.double(dya), 
            nb=as.integer(nb),
            x0b=as.double(x0b),
            y0b=as.double(y0b),
            dxb=as.double(dxb),
            dyb=as.double(dyb), 
            eps=as.double(eps),
            ok=as.integer(integer(1)))
  hit <- (out$ok != 0)
  return(hit)
}

selfcrossing.psp <- function(A) {
  verifyclass(A, "psp")
  eps <- .Machine$double.eps

  n <- A$n
  eA <- A$ends
  x0 <- eA$x0
  y0 <- eA$y0
  dx <- eA$x1 - eA$x0
  dy <- eA$y1 - eA$y0

  useCall <- spatstat.options("selfcrossing.psp.useCall")
  if(!useCall) {
    # old C routine
    out <- .C("xysegXint",
              n=as.integer(n),
              x0=as.double(x0),
              y0=as.double(y0),
              dx=as.double(dx),
              dy=as.double(dy), 
              eps=as.double(eps),
              xx=as.double(numeric(n^2)),
              yy=as.double(numeric(n^2)),
              ti=as.double(numeric(n^2)),
              tj=as.double(numeric(n^2)),
              ok=as.integer(integer(n^2)))

    ok <- (matrix(out$ok, n, n) != 0)
    xx <- matrix(out$xx, n, n)
    yy <- matrix(out$yy, n, n)
    xx <- as.vector(xx[ok])
    yy <- as.vector(yy[ok])
  } else {
    # new
    storage.mode(x0) <- storage.mode(y0) <- "double"
    storage.mode(dx) <- storage.mode(dy) <- "double"
    storage.mode(eps) <- "double"
    out <- .Call("CxysegXint",
                 x0, 
                 y0, 
                 dx, 
                 dy, 
    	         eps)
#                 PACKAGE="spatstat")
    xx <- out[[5]]
    yy <- out[[6]]
  }
  result <- ppp(xx, yy, window=A$window, check=FALSE)
  return(result)
}


test.selfcrossing.psp <- function(A) {
  verifyclass(A, "psp")
  eps <- .Machine$double.eps

  n <- A$n
  eA <- A$ends
  x0 <- eA$x0
  y0 <- eA$y0
  dx <- eA$x1 - eA$x0
  dy <- eA$y1 - eA$y0

  out <- .C("xysxi",
            na=as.integer(n),
            x0=as.double(x0),
            y0=as.double(y0),
            dx=as.double(dx),
            dy=as.double(dy), 
            eps=as.double(eps),
            ok=as.integer(integer(n*n)))
  hit <- (matrix(out$ok, n, n) != 0)
  return(hit)
}

selfcut.psp <- function(A, ..., eps) {
  stopifnot(is.psp(A))
#  n <- A$n
  eA <- A$ends
  x0 <- eA$x0
  y0 <- eA$y0
  dx <- eA$x1 - eA$x0
  dy <- eA$y1 - eA$y0
  if(missing(eps) || is.null(eps)) {
    eps <- sqrt(.Machine$double.eps) * diameter(Frame(A))
  } else {
    check.1.real(eps)
    stopifnot(eps >= 0)
  }
  ## identify self-crossings
  eps <- .Machine$double.eps
  storage.mode(x0) <- storage.mode(y0) <- "double"
  storage.mode(dx) <- storage.mode(dy) <- "double"
  storage.mode(eps) <- "double"
  zz <- .Call("CxysegXint",
              x0, 
              y0, 
              dx, 
              dy, 
              eps)
  if(length(zz[[1]]) == 0)
    return(A)
  ##
  names(zz) <- c("i", "j", "ti", "tj", "x", "y")
  df <- as.data.frame(zz)
  df$i <- df$i + 1L
  df$j <- df$j + 1L
  ##
  gone <- with(df, unique(c(i,j)))
  newends <- as.matrix(eA)
  newends <- newends[-gone, , drop=FALSE]
  newmarx <- marx <- marks(A)
  if(mama <- !is.null(marx)) 
    newmarx <- as.data.frame(marx)[-gone, ,drop=FALSE]
  ## cut each segment using the *provided* values of x,y
  for(ii in gone) {
    ## assemble cuts through segment ii
    imatch <- with(df, which(i == ii))
    jmatch <- with(df, which(j == ii))
    df.i <- with(df,
                 data.frame(t=c(ti[imatch], tj[jmatch]),
                            x=x[c(imatch, jmatch)],
                            y=y[c(imatch, jmatch)]))
    # discard T-junctions
    ok <- with(df.i, t > 0 & t < 1)
    df.i <- df.i[ok, ,drop=FALSE]
    # order the pieces
    ord <- with(df.i, order(t))
    df.i <- df.i[ord, , drop=FALSE]
    ## add endpoints
    xnew <- c(eA[ii,"x0"], df.i$x, eA[ii,"x1"])
    ynew <- c(eA[ii,"y0"], df.i$y, eA[ii,"y1"])
    m <- length(xnew)
    newsegs <- cbind(xnew[-m], ynew[-m], xnew[-1], ynew[-1])
    newends <- rbind(newends, newsegs)
    if(mama)
      newmarx <- rbind(newmarx, marx[rep(ii, m-1), , drop=FALSE])
  }
  Y <- as.psp(newends, window=Window(A), marks=newmarx)
  if(eps > 0) {
    ok <- (lengths.psp(Y) > eps)
    if(any(!ok)) Y <- Y[ok]
  }
  return(Y)
}
  
