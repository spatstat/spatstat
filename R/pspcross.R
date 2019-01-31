#
#    pspcross.R
#
#    Intersections of line segments
#    
#    $Revision: 1.24 $   $Date: 2019/01/31 08:11:34 $
#
#
crossing.psp <- function(A,B,fatal=TRUE,details=FALSE) {
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
              ok=as.integer(integer(na * nb)),
              PACKAGE = "spatstat")
    
    ok <- (matrix(out$ok, na, nb) != 0)
    xx <- matrix(out$xx, na, nb)
    yy <- matrix(out$yy, na, nb)
    xx <- as.vector(xx[ok])
    yy <- as.vector(yy[ok])
    if(details) {
      ia <- as.vector(row(ok)[ok])
      jb <- as.vector(col(ok)[ok])
      ta <- as.vector(matrix(out$ta, na, nb)[ok])
      tb <- as.vector(matrix(out$tb, na, nb)[ok])
    }
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
    	           eps,
                 PACKAGE="spatstat")
    xx <- out[[5]]
    yy <- out[[6]]
    if(details) {
      ia <- out[[1L]] + 1L
      jb <- out[[2L]] + 1L
      ta <- out[[3L]]
      tb <- out[[4L]]
    }
  }
  result <- ppp(xx, yy, window=ABW, check=FALSE)
  if(details)
    marks(result) <- data.frame(iA=ia, jB=jb, tA=ta, tB=tb)
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
            ok=as.integer(integer(na * nb)),
            PACKAGE = "spatstat")

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
            ok=as.integer(integer(1L)),
            PACKAGE = "spatstat")
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
              ok=as.integer(integer(n^2)),
              PACKAGE = "spatstat")

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
    	           eps,
                 PACKAGE="spatstat")
    xx <- out[[5L]]
    yy <- out[[6L]]
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
            ok=as.integer(integer(n*n)),
            PACKAGE = "spatstat")
  hit <- (matrix(out$ok, n, n) != 0)
  return(hit)
}

selfcut.psp <- function(A, ..., eps) {
  stopifnot(is.psp(A))
  n <- A$n
  eA <- A$ends
  x0 <- eA$x0
  y0 <- eA$y0
  x1 <- eA$x1
  y1 <- eA$y1
  dx <- x1 - x0
  dy <- y1 - y0
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
              eps,
              PACKAGE = "spatstat")
  if(length(zz[[1]]) == 0) {
    ## no dissection required
    attr(A, "camefrom") <- seq_len(n)
    return(A)
  }
  ##
  names(zz) <- c("i", "j", "ti", "tj", "x", "y")
  icross  <- zz$i + 1L
  jcross  <- zz$j + 1L
  ticross <- zz$ti
  tjcross <- zz$tj
  xcross  <- zz$x
  ycross  <- zz$y
  ## which segments are split...
  gone <- unique(c(icross, jcross))
  ## ... and which are not
  retained <- setdiff(seq_len(n), gone)
  ## initialise result
  ## start with all segments which are retained
  x0out <- x0[retained]
  y0out <- y0[retained]
  x1out <- x1[retained]
  y1out <- y1[retained]
  camefrom <- retained
  ## cut each segment using the *provided* values of x,y
  for(ii in gone) {
    ## assemble cuts through segment ii
    imatch <- which(icross == ii)
    jmatch <- which(jcross == ii)
    ijmatch <- c(imatch, jmatch)
    tt <- c(ticross[imatch], tjcross[jmatch])
    xx <- xcross[ijmatch]
    yy <- ycross[ijmatch]
    # discard T-junctions
    ok <- (tt > 0 & tt < 1)
    tt <- tt[ok]
    xx <- xx[ok]
    yy <- yy[ok]
    # order the pieces
    ord <- order(tt)
    xx <- xx[ord]
    yy <- yy[ord]
    ## add endpoints of old segment
    xnew <- c(x0[ii], xx, x1[ii])
    ynew <- c(y0[ii], yy, y1[ii])
    ## append to result
    m <- length(xnew)
    x0out <- c(x0out, xnew[-m])
    y0out <- c(y0out, ynew[-m])
    x1out <- c(x1out, xnew[-1L])
    y1out <- c(y1out, ynew[-1L])
    camefrom <- c(camefrom, rep(ii, m-1L))
  }
  marx <- marks(A)
  marxout <- if(is.null(marx)) NULL else
             as.data.frame(marx)[camefrom, , drop=FALSE]
  Y <- psp(x0out, y0out, x1out, y1out, window=Window(A), marks=marxout)
  if(eps > 0) {
    ok <- (lengths.psp(Y) > eps)
    if(!all(ok)) {
      Y <- Y[ok]
      camefrom <- camefrom[ok]
    }
  }
  attr(Y, "camefrom") <- camefrom
  return(Y)
}
  
