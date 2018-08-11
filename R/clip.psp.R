#
# clip.psp.R
#
#    $Revision: 1.23 $   $Date: 2018/01/23 02:40:14 $
#
#
 
########################################################
# clipping operation (for subset)
########################################################

clip.psp <- function(x, window, check=TRUE, fragments=TRUE) {
  verifyclass(x, "psp")
  verifyclass(window, "owin")
  if(is.vanilla(unitname(window)))
     unitname(window) <- unitname(x)
  if(check && !is.subset.owin(window, x$window))
    warning("The clipping window is not a subset of the window containing the line segment pattern x")
  if(x$n == 0) {
    emptypattern <- psp(numeric(0), numeric(0), numeric(0), numeric(0),
                      window=window, marks=x$marks)
    return(emptypattern)
  }
  switch(window$type,
         rectangle={
           result <- cliprect.psp(x, window, fragments=fragments)
         },
         polygonal={
           result <- clippoly.psp(x, window, fragments=fragments)
         },
         mask={
           result <- clippoly.psp(x, as.polygonal(window), fragments=fragments)
           result$window <- window
         })
  return(result)
}


#####
#
#  clipping to a rectangle
#
cliprect.psp <- local({
  
  cliprect.psp <- function(x, window, fragments=TRUE) {
    verifyclass(x, "psp")
    verifyclass(window, "owin")
    ends <- x$ends
    marx <- marks(x, dfok=TRUE)
    #' find segments which are entirely inside the window
    #' (by convexity)
    in0 <- inside.owin(ends$x0, ends$y0, window)
    in1 <- inside.owin(ends$x1, ends$y1, window)
    ok <- in0 & in1
    #' if all segments are inside, return them
    if(all(ok))
      return(as.psp(ends, window=window, marks=marx, check=FALSE))
    #' otherwise, store those segments which are inside the window
    ends.inside <- ends[ok, , drop=FALSE]
    marks.inside <- marx %msub% ok
    x.inside <- as.psp(ends.inside, window=window, marks=marks.inside,
                       check=FALSE)
    if(!fragments)
      return(x.inside)
    #' now consider the rest
    ends <- ends[!ok, , drop=FALSE]
    in0 <- in0[!ok] 
    in1 <- in1[!ok]
    marx <- marx %msub% (!ok)
    #' first clip segments to the range x \in [xmin, xmax]
    #' use parametric coordinates
    tx <- cbind(ifelse0NA(between(ends$x0, window$xrange)),
                ifelse1NA(between(ends$x1, window$xrange)),
                tvalue(ends$x0, ends$x1, window$xrange[1L]),
                tvalue(ends$x0, ends$x1, window$xrange[2L]))
    #' discard segments which do not lie in the x range 
    nx <- apply(!is.na(tx), 1L, sum)
    ok <- (nx >= 2)
    if(!any(ok))
      return(x.inside)
    ends <- ends[ok, , drop=FALSE]
    tx   <- tx[ok, , drop=FALSE]
    in0  <- in0[ok]
    in1  <- in1[ok]
    marx <- marx %msub% ok
    #' Clip the segments to the x range
    tmin <- apply(tx, 1L, min, na.rm=TRUE)
    tmax <- apply(tx, 1L, max, na.rm=TRUE)
    dx <- ends$x1 - ends$x0
    dy <- ends$y1 - ends$y0
    ends.xclipped <- data.frame(x0=ends$x0 + tmin * dx,
                                y0=ends$y0 + tmin * dy,
                                x1=ends$x0 + tmax * dx,
                                y1=ends$y0 + tmax * dy)
    #' Now clip the segments to the range y \in [ymin, ymax]
    ends <- ends.xclipped
    in0 <- inside.owin(ends$x0, ends$y0, window)
    in1 <- inside.owin(ends$x1, ends$y1, window)
    ty <- cbind(ifelse0NA(in0),
                ifelse1NA(in1),
                tvalue(ends$y0, ends$y1, window$yrange[1L]),
                tvalue(ends$y0, ends$y1, window$yrange[2L]))
    #' discard segments which do not lie in the y range 
    ny <- apply(!is.na(ty), 1L, sum)
    ok <- (ny >= 2)
    if(!any(ok))
      return(x.inside)
    ends <- ends[ok, , drop=FALSE]
    ty   <- ty[ok, , drop=FALSE]
    in0  <- in0[ok]
    in1  <- in1[ok]
    marx <- marx %msub% ok
    #' Clip the segments to the y range
    tmin <- apply(ty, 1L, min, na.rm=TRUE)
    tmax <- apply(ty, 1L, max, na.rm=TRUE)
    dx <- ends$x1 - ends$x0
    dy <- ends$y1 - ends$y0
    ends.clipped <- data.frame(x0=ends$x0 + tmin * dx,
                               y0=ends$y0 + tmin * dy,
                               x1=ends$x0 + tmax * dx,
                               y1=ends$y0 + tmax * dy)
    marks.clipped <- marx
    #' OK - segments clipped
    #' Put them together with the unclipped ones
    ends.all <- rbind(ends.inside, ends.clipped)
    marks.all <- marks.inside %mapp% marks.clipped
    as.psp(ends.all, window=window, marks=marks.all, check=FALSE)
  }

  small <- function(x) { abs(x) <= .Machine$double.eps }

  tvalue <- function(z0, z1, zt) {
    y1 <- z1 - z0
    yt <- zt - z0
    tval <- ifelseAX(small(y1), 0.5, yt/y1)
    betwee <- (yt * (zt - z1)) <= 0
    result <- ifelseXB(betwee, tval, NA)
    return(result)
  }

  between <- function(x, r) { ((x-r[1L]) * (x-r[2L])) <= 0 }

  cliprect.psp
})

############################
#
#   clipping to a polygonal window
#

clippoly.psp <- function(s, window, fragments=TRUE) {
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
  x1s <- es$x1
  y1s <- es$y1
  dxs <- x1s - x0s
  dys <- y1s - y0s

  bdry <- edges(window)
  nw <- bdry$n
  ew <- bdry$ends
  x0w <- ew$x0
  y0w <- ew$y0
  dxw <- ew$x1 - ew$x0
  dyw <- ew$y1 - ew$y0

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
            PACKAGE = "spatstat")

  hitting <- (matrix(out$ok, ns, nw) != 0)
  ts <- matrix(out$ta, ns, nw)

  anyhit <- matrowany(hitting)
  
  if(!fragments) {
    #' retain only segments that avoid the boundary entirely
    leftin <- inside.owin(es$x0, es$y0, window)
    rightin <- inside.owin(es$x1, es$y1, window)
    ok <- !anyhit & leftin & rightin
    return(as.psp(es[ok,,drop=FALSE],
                  window=window,
		  marks=marx %msub% ok,
		  check=FALSE))
  }
  # form all the chopped segments (whether in or out)

  #' initially empty
  chopx0 <- chopy0 <- chopx1 <- chopy1 <- numeric(0)
  chopmarks <- marx %msub% integer(0)
    
  for(seg in seq_len(ns)) {
    #' coordinates of segment number 'seg'
    segx0 <- x0s[seg]
    segy0 <- y0s[seg]
    segx1 <- x1s[seg]
    segy1 <- y1s[seg]
    if(has.marks)
      segmarks <- marx %msub% seg
    if(!anyhit[seg]) {
      #' no intersection with boundary - add single segment
      chopx0 <- c(chopx0, segx0)
      chopy0 <- c(chopy0, segy0)
      chopx1 <- c(chopx1, segx1)
      chopy1 <- c(chopy1, segy1)
      if(has.marks)
        chopmarks <- chopmarks %mapp% segmarks
    } else {
      #' crosses boundary - add several pieces
      tvals <- ts[seg,]
      tvals <- sort(tvals[hitting[seg,]])
      dx <- segx1 - segx0
      dy <- segy1 - segy0
      chopx0 <- c(chopx0, segx0 + c(0,tvals) * dx)
      chopy0 <- c(chopy0, segy0 + c(0,tvals) * dy)
      chopx1 <- c(chopx1, segx0 + c(tvals,1) * dx)
      chopy1 <- c(chopy1, segy0 + c(tvals,1) * dy)
      if(has.marks) {
        npieces <- length(tvals) + 1L
        chopmarks <- chopmarks %mapp% (segmarks %mrep% npieces)
      }
    }
  }

  chopped <- psp(chopx0, chopy0, chopx1, chopy1,
                 window=boundingbox(Window(s), window),
                 marks=chopmarks)
  
  # select those chopped segments which are inside the window
  mid <- midpoints.psp(chopped)
  ins <- inside.owin(mid$x, mid$y, window)
  retained <- chopped[ins]
  retained$window <- window
  return(retained)
}


