#
# clip.psp.R
#
#    $Revision: 1.15 $   $Date: 2013/05/01 05:45:01 $
#
#
 
########################################################
# clipping operation (for subset)
########################################################

clip.psp <- function(x, window, check=TRUE) {
  verifyclass(x, "psp")
  verifyclass(window, "owin")
  if(check && !is.subset.owin(window, x$window))
    warning("The clipping window is not a subset of the window containing the line segment pattern x")
  if(x$n == 0) {
    emptypattern <- psp(numeric(0), numeric(0), numeric(0), numeric(0),
                      window=window, marks=x$marks)
    return(emptypattern)
  }
  switch(window$type,
         rectangle=cliprect.psp(x, window),
         polygonal=clippoly.psp(x, window),
         mask=stop("sorry, clipping is not implemented for masks"))
}


#####
#
#  clipping to a rectangle
#
cliprect.psp <- function(x, window) {
  verifyclass(x, "psp")
  verifyclass(window, "owin")
  ends <- x$ends
  marx <- marks(x, dfok=TRUE)
  # find segments which are entirely inside the window
  # (by convexity)
  in0 <- inside.owin(ends$x0, ends$y0, window)
  in1 <- inside.owin(ends$x1, ends$y1, window)
  ok <- in0 & in1
  # if all segments are inside, return them
  if(all(ok))
    return(as.psp(ends, window=window, marks=marx, check=FALSE))
  # otherwise, store those segments which are inside the window
  ends.inside <- ends[ok, , drop=FALSE]
  marks.inside <- marx %msub% ok
  x.inside <- as.psp(ends.inside, window=window, marks=marks.inside, check=FALSE)
  # now consider the rest
  ends <- ends[!ok, , drop=FALSE]
  in0 <- in0[!ok] 
  in1 <- in1[!ok]
  marx <- marx %msub% (!ok)
  # first clip segments to the range x \in [xmin, xmax]
  # use parametric coordinates
  small <- function(x) { abs(x) <= .Machine$double.eps }
  tvalue <- function(z0, z1, zt) {
    y1 <- z1 - z0
    yt <- zt - z0
    tval <- ifelseAX(small(y1), 0.5, yt/y1)
    betwee <- (yt * (zt - z1)) <= 0
    result <- ifelseXB(betwee, tval, NA)
    return(result)
  }
  between <- function(x, r) { ((x-r[1]) * (x-r[2])) <= 0 }
  tx <- cbind(ifelse0NA(between(ends$x0, window$xrange)),
              ifelse1NA(between(ends$x1, window$xrange)),
              tvalue(ends$x0, ends$x1, window$xrange[1]),
              tvalue(ends$x0, ends$x1, window$xrange[2]))
  # discard segments which do not lie in the x range 
  nx <- apply(!is.na(tx), 1, sum)
  ok <- (nx >= 2)
  if(!any(ok))
    return(x.inside)
  ends <- ends[ok, , drop=FALSE]
  tx   <- tx[ok, , drop=FALSE]
  in0  <- in0[ok]
  in1  <- in1[ok]
  marx <- marx %msub% ok
  # Clip the segments to the x range
  tmin <- apply(tx, 1, min, na.rm=TRUE)
  tmax <- apply(tx, 1, max, na.rm=TRUE)
  dx <- ends$x1 - ends$x0
  dy <- ends$y1 - ends$y0
  ends.xclipped <- data.frame(x0=ends$x0 + tmin * dx,
                             y0=ends$y0 + tmin * dy,
                             x1=ends$x0 + tmax * dx,
                             y1=ends$y0 + tmax * dy)
  # Now clip the segments to the range y \in [ymin, ymax]
  ends <- ends.xclipped
  in0 <- inside.owin(ends$x0, ends$y0, window)
  in1 <- inside.owin(ends$x1, ends$y1, window)
  ty <- cbind(ifelse0NA(in0),
              ifelse1NA(in1),
              tvalue(ends$y0, ends$y1, window$yrange[1]),
              tvalue(ends$y0, ends$y1, window$yrange[2]))
  # discard segments which do not lie in the y range 
  ny <- apply(!is.na(ty), 1, sum)
  ok <- (ny >= 2)
  if(!any(ok))
    return(x.inside)
  ends <- ends[ok, , drop=FALSE]
  ty   <- ty[ok, , drop=FALSE]
  in0  <- in0[ok]
  in1  <- in1[ok]
  marx <- marx %msub% ok
  # Clip the segments to the y range
  tmin <- apply(ty, 1, min, na.rm=TRUE)
  tmax <- apply(ty, 1, max, na.rm=TRUE)
  dx <- ends$x1 - ends$x0
  dy <- ends$y1 - ends$y0
  ends.clipped <- data.frame(x0=ends$x0 + tmin * dx,
                             y0=ends$y0 + tmin * dy,
                             x1=ends$x0 + tmax * dx,
                             y1=ends$y0 + tmax * dy)
  marks.clipped <- marx
  # OK - segments clipped
  # Put them together with the unclipped ones
  ends.all <- rbind(ends.inside, ends.clipped)
  marks.all <- marks.inside %mapp% marks.clipped
  as.psp(ends.all, window=window, marks=marks.all, check=FALSE)
}


############################
#
#   clipping to a polygonal window
#

clippoly.psp <- function(s, window) {
  verifyclass(s, "psp")
  verifyclass(window, "owin")
  stopifnot(window$type == "polygonal")
  marx <- marks(s)
  has.marks <- !is.null(marx)
  
  eps <- .Machine$double.eps

  # find the intersection points between segments and window edges
  
  ns <- s$n
  es <- s$ends
  x0s <- es$x0
  y0s <- es$y0
  dxs <- es$x1 - es$x0
  dys <- es$y1 - es$y0

  bdry <- as.psp(window)
  nw <- bdry$n
  ew <- bdry$ends
  x0w <- ew$x0
  y0w <- ew$y0
  dxw <- ew$x1 - ew$x0
  dyw <- ew$y1 - ew$y0

  DUP <- spatstat.options("dupC")
  out <- .C("xysegint",
            na=as.integer(ns),
            x0a=as.double(x0s),
            y0a=as.double(y0s),
            dxa=as.double(dxs),
            dya=as.double(dys), 
            nb=as.integer(nw),
            x0b=as.double(x0w),
            y0b=as.double(y0w),
            dxb=as.double(dxw),
            dyb=as.double(dyw), 
            eps=as.double(eps),
            xx=as.double(numeric(ns * nw)),
            yy=as.double(numeric(ns * nw)),
            ta=as.double(numeric(ns * nw)),
            tb=as.double(numeric(ns * nw)),
            ok=as.integer(integer(ns * nw)),
     DUP=DUP)
#     PACKAGE="spatstat")

  ok <- (matrix(out$ok, ns, nw) != 0)
  ts <- matrix(out$ta, ns, nw)

  # form all the chopped segments (whether in or out)

  chopped <- empty <- s[numeric(0)]
  chopped$window <- bounding.box(s$window, window)
    
  for(seg in seq_len(ns)) {
    segment <- s$ends[seg, , drop=FALSE]
    hit <- ok[seg, ]
    if(!any(hit)) {
      # no intersection with boundary - add single segment
      chopped$ends <- rbind(chopped$ends, segment)
    if(has.marks) chopped$marks <- (chopped$marks) %mapp% (marx %msub% seg)
    } else {
      # crosses boundary - add several pieces
      tvals <- ts[seg,]
      tvals <- sort(tvals[hit])
      x0 <- segment$x0
      dx <- segment$x1 - x0
      y0 <- segment$y0
      dy <- segment$y1 - y0
      newones <- data.frame(x0 = x0 + c(0,tvals) * dx,
                            y0 = y0 + c(0,tvals) * dy,
                            x1 = x0 + c(tvals,1) * dx,
                            y1 = y0 + c(tvals,1) * dy)
      chopped$ends <- rbind(chopped$ends, newones)
      if(has.marks) {
        hitmarks <- marx %msub% seg
        newmarks <- hitmarks %mrep% nrow(newones)
        chopped$marks <-  (chopped$marks) %mapp% newmarks
      }
    }
  }
  chopped$n <- nrow(chopped$ends)
  
  # select those chopped segments which are inside the window
  mid <- midpoints.psp(chopped)
  ins <- inside.owin(mid$x, mid$y, window)
  retained <- chopped[ins]
  retained$window <- window
  return(retained)
}


