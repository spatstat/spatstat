#
#
#   disc.R
#
#   $Revision: 1.28 $ $Date: 2018/03/23 01:26:43 $
#
#   Compute the disc of radius r in a linear network
#
#   
lineardisc <- function(L, x=locator(1), r, plotit=TRUE,
                       cols=c("blue", "red", "green")) {
  # L is the linear network (object of class "linnet")
  # x is the centre point of the disc
  # r is the radius of the disc
  #
  stopifnot(inherits(L, "linnet"))
  check.1.real(r)
  if(L$sparse) {
    message("Converting linear network to non-sparse representation..")
    L <- as.linnet(L, sparse=FALSE)
  }
  lines <- L$lines
  vertices <- L$vertices
  lengths <- lengths.psp(lines)
  win <- L$window
  #
  # project x to nearest segment
  if(missing(x))
    x <- clickppp(1, win, add=TRUE)
  else
    x <- as.ppp(x, win)
  pro <- project2segment(x, lines)
  # which segment?
  startsegment <- pro$mapXY
  # parametric position of x along this segment
  startfraction <- pro$tp
  # vertices at each end of this segment
  A <- L$from[startsegment]
  B <- L$to[startsegment]
  # distances from x to  A and B
  dxA <- startfraction * lengths[startsegment]
  dxB <- (1-startfraction) * lengths[startsegment]
  # is r large enough to reach both A and B?
  startfilled <- (max(dxA, dxB) <= r)
  # compute vector of shortest path distances from x to each vertex j,
  # going through A:
  dxAv <- dxA + L$dpath[A,]
  # going through B:
  dxBv <- dxB + L$dpath[B,]
  # going either through A or through B:
  dxv <- pmin.int(dxAv, dxBv)
  # Thus dxv[j] is the shortest path distance from x to vertex j.
  #
  # Determine which vertices are inside the disc of radius r
  covered <- (dxv <= r)
  # Thus covered[j] is TRUE if the j-th vertex is inside the disc.
  #
  # Determine which line segments are completely inside the disc
  #
  from <- L$from
  to   <- L$to
  # ( a line segment is inside the disc if the shortest distance
  #   from x to one of its endpoints, plus the length of the segment,
  #   is less than r ....
  allinside <- (dxv[from] + lengths <= r) | (dxv[to] + lengths <= r)
  #   ... or alternatively, if the sum of the
  #   two residual distances exceeds the length of the segment )
  residfrom <- pmax.int(0, r - dxv[from])
  residto   <- pmax.int(0, r - dxv[to])
  allinside <- allinside | (residfrom + residto >= lengths)
  # start segment is special
  allinside[startsegment] <- startfilled
  # Thus allinside[k] is TRUE if the k-th segment is inside the disc
  
  # Collect all these segments
  disclines <- lines[allinside]
  #
  # Determine which line segments cross the boundary of the disc
  boundary <- (covered[from] | covered[to]) & !allinside
  # For each of these, calculate the remaining distance at each end
  resid.from <- ifelseXB(boundary, pmax.int(r - dxv[from], 0), 0)
  resid.to   <- ifelseXB(boundary, pmax.int(r - dxv[to],   0), 0)
  # Where the remaining distance is nonzero, create segment and endpoint
  okfrom <- (resid.from > 0)
  okfrom[startsegment] <- FALSE
  if(any(okfrom)) {
    v0 <- vertices[from[okfrom]]
    v1 <- vertices[to[okfrom]]
    tp <- (resid.from/lengths)[okfrom]
    vfrom <- ppp((1-tp)*v0$x + tp*v1$x,
                 (1-tp)*v0$y + tp*v1$y,
                 window=win)
    extralinesfrom <- as.psp(from=v0, to=vfrom)
  } else vfrom <- extralinesfrom <- NULL
  #
  okto <- (resid.to > 0)
  okto[startsegment] <- FALSE
  if(any(okto)) {
    v0 <- vertices[to[okto]]
    v1 <- vertices[from[okto]]
    tp <- (resid.to/lengths)[okto]
    vto <- ppp((1-tp)*v0$x + tp*v1$x,
               (1-tp)*v0$y + tp*v1$y,
               window=win)
    extralinesto <- as.psp(from=v0, to=vto)
  } else vto <- extralinesto <- NULL
  #
  # deal with special case where start segment is not fully covered
  if(!startfilled) {
    vA <- vertices[A]
    vB <- vertices[B]
    rfrac <- r/lengths[startsegment]
    tleft <- pmax.int(startfraction-rfrac, 0)
    tright <- pmin.int(startfraction+rfrac, 1)
    vleft <- ppp((1-tleft) * vA$x + tleft * vB$x,
                 (1-tleft) * vA$y + tleft * vB$y,
                 window=win)
    vright <- ppp((1-tright) * vA$x + tright * vB$x,
                  (1-tright) * vA$y + tright * vB$y,
                  window=win)
    startline <- as.psp(from=vleft, to=vright)
    startends <- superimpose(if(!covered[A]) vleft else NULL,
                             if(!covered[B]) vright else NULL)
  } else startline <- startends <- NULL
  #
  # combine all lines
  disclines <- superimpose(disclines,
                           extralinesfrom, extralinesto, startline,
                           W=win, check=FALSE)
  # combine all disc endpoints
  discends <- superimpose(vfrom, vto, vertices[dxv == r], startends,
                          W=win, check=FALSE)
  #
  if(plotit) {
    if(dev.cur() == 1) {
      # null device - initialise a plot
      plot(L, main="")
    }
    points(x, col=cols[1L], pch=16)
    plot(disclines, add=TRUE, col=cols[2L], lwd=2)
    plot(discends, add=TRUE, col=cols[3L], pch=16)
  }
  return(list(lines=disclines, endpoints=discends))
}

countends <- function(L, x=locator(1), r, toler=NULL) {
  # L is the linear network (object of class "linnet")
  # x is the centre point of the disc
  # r is the radius of the disc
  #
  stopifnot(inherits(L, "linnet"))
  sparse <- L$sparse %orifnull% is.null(L$dpath)
  if(sparse)
    stop(paste("countends() does not support linear networks",
               "that are stored in sparse matrix format.",
               "Please convert the data using as.linnet(sparse=FALSE)"),
         call.=FALSE)
  # get x
  if(missing(x))
    x <- clickppp(1, Window(L), add=TRUE)
  if(!inherits(x, "lpp"))
    x <- as.lpp(x, L=L)
  np <- npoints(x)
  
  if(length(r) != np)
    stop("Length of vector r does not match number of points in x")
  #
  if(!is.connected(L)) {
    #' disconnected network - split into components
    result <- numeric(np)
    lab <- connected(L, what="labels")
    subsets <- split(seq_len(nvertices(L)), factor(lab))
    for(subi in subsets) {
      xi <- thinNetwork(x, retainvertices=subi)
      witch <- which(attr(xi, "retainpoints"))
      ok <- is.finite(r[witch])
      witchok <- witch[ok]
      result[witchok] <-
        countends(domain(xi), xi[ok], r[witchok], toler=toler)      
    }
    return(result)
  }
  lines <- L$lines
  vertices <- L$vertices
  lengths <- lengths.psp(lines)
  dpath <- L$dpath
  nv <- vertices$n
  ns <- lines$n
  #
  if(!spatstat.options("Ccountends")) {
    #' interpreted code
    result <- integer(np)
    for(i in seq_len(np)) 
      result[i] <- npoints(lineardisc(L, x[i], r[i], plotit=FALSE)$endpoints)
    return(result)
  }
  # extract coordinates
  coo <- coords(x)
  #' which segment
  startsegment <- coo$seg 
  # parametric position of x along this segment
  startfraction <- coo$tp
  # convert indices to C 
  seg0 <- startsegment - 1L
  from0 <- L$from - 1L
  to0   <- L$to - 1L
  # determine numerical tolerance
  if(is.null(toler)) {
    toler <- default.linnet.tolerance(L)
  } else {
    check.1.real(toler)
    stopifnot(toler > 0)
  }
  zz <- .C("Ccountends",
           np = as.integer(np),
           f = as.double(startfraction),
           seg = as.integer(seg0),
           r = as.double(r), 
           nv = as.integer(nv), 
           xv = as.double(vertices$x),
           yv = as.double(vertices$y),  
           ns = as.integer(ns),
           from = as.integer(from0),
           to = as.integer(to0), 
           dpath = as.double(dpath),
           lengths = as.double(lengths),
           toler=as.double(toler),
           nendpoints = as.integer(integer(np)),
           PACKAGE = "spatstat")
  zz$nendpoints
}

default.linnet.tolerance <- function(L) {
  # L could be a linnet or psp
  if(!is.null(toler <- L$toler)) return(toler)
  len2 <- lengths.psp(as.psp(L), squared=TRUE)
  len2pos <- len2[len2 > 0]
  toler <- if(length(len2pos) == 0) 0 else (0.001 * sqrt(min(len2pos)))
  toler <- makeLinnetTolerance(toler)
  return(toler)
}

makeLinnetTolerance <- function(toler) {
  max(sqrt(.Machine$double.xmin),
      toler[is.finite(toler)], na.rm=TRUE)
}




  
