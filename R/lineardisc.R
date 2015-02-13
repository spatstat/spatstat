#
#
#   disc.R
#
#   $Revision: 1.19 $ $Date: 2014/11/10 05:39:57 $
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
    points(x, col=cols[1], pch=16)
    plot(disclines, add=TRUE, col=cols[2], lwd=2)
    plot(discends, add=TRUE, col=cols[3], pch=16)
  }
  return(list(lines=disclines, endpoints=discends))
}

countends <- function(L, x=locator(1), r) {
  # L is the linear network (object of class "linnet")
  # x is the centre point of the disc
  # r is the radius of the disc
  #
  stopifnot(inherits(L, "linnet"))
  lines <- L$lines
  vertices <- L$vertices
  lengths <- lengths.psp(lines)
  dpath <- L$dpath
  win <- L$window
  nv <- vertices$n
  ns <- lines$n
  # get x
  if(missing(x))
    x <- clickppp(1, win, add=TRUE)
  else
    x <- as.ppp(x, win)
  #
  np <- npoints(x)
  if(length(r) != np)
    stop("Length of vector r does not match number of points in x")
  # project x to nearest segment
  pro <- project2segment(x, lines)
  # which segment?
  startsegment <- pro$mapXY
  # parametric position of x along this segment
  startfraction <- pro$tp

  # convert indices to C 
  seg0 <- startsegment - 1L
  from0 <- L$from - 1L
  to0   <- L$to - 1L
  toler <- 0.001 * min(lengths)
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
           nendpoints = as.integer(integer(np)))
  zz$nendpoints
}
