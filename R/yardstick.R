##
##   yardstick.R
##
##   Simple class 'yardstick' to display scale information
##
##   $Revision: 1.2 $ $Date: 2014/11/20 03:21:07 $

yardstick <- function(x0, y0, x1, y1, txt) {
  nomore <- missing(y0) && missing(x1) && missing(y1) 
  if(is.ppp(x0) && nomore) {
    if(npoints(x0) != 2) stop("x0 should consist of exactly 2 points")
    X <- x0
  } else if(is.psp(x0) && nomore) {
    if(nobjects(x0) != 1) stop("x0 should consist of exactly 1 segment")
    X <- endpoints.psp(x0)
  } else {
    xx <- c(x0, x1)
    yy <- c(y0, y1)
    B <- boundingbox(list(x=xx, y=yy))
    X <- ppp(xx, yy, window=grow.rectangle(B, diameter(B)))
  }
  attr(X, "txt") <- txt
  class(X) <- c("yardstick", class(X))
  return(X)
}

"[.yardstick" <- function(x, ...) {
  yardstick(NextMethod("["), txt=attr(x, "txt"))
}

shift.yardstick <- function(X, ...) {
  yardstick(NextMethod("shift"), txt=attr(X, "txt"))
}

plot.yardstick <- local({

  myarrows <- function(x0, y0, x1, y1, ...,
                       left=TRUE, right=TRUE,
                       angle=20, frac=0.25,
                       main, show.all, add) {
    segments(x0, y0, x1, y1, ...)
    if(left || right) {
      ang <- angle * pi/180
      co <- cos(ang)
      si <- sin(ang)
      dx <- x1-x0
      dy <- y1-y0
      le <- sqrt(dx^2 + dy^2)
      rot <- matrix(c(dx, dy, -dy, dx)/le, 2, 2)
      arlen <- frac * le
      up <- arlen * (rot %*% c(co, si))
      lo <- arlen * (rot %*% c(co, -si))
      if(left) {
        segments(x0, y0, x0+up[1], y0+up[2], ...)
        segments(x0, y0, x0+lo[1], y0+lo[2], ...)
      }
      if(right) {
        segments(x1, y1, x1-up[1], y1-up[2], ...)
        segments(x1, y1, x1-lo[1], y1-lo[2], ...)
      }
    }
    return(invisible(NULL))
  }

  plot.yardstick <- function(x, ...,
                             angle=20,
                             frac=1/8,
                             cex=1,
                             pos=NULL,
                             split=FALSE,
                             shrink=1/4,
                             do.plot=TRUE) {
    if(do.plot) {
      A <- as.numeric(coords(x)[1,])
      B <- as.numeric(coords(x)[2,])
      M <- (A+B)/2
      if(!split) {
        ## double-headed arrow
        myarrows(A[1], A[2], B[1], B[2], ..., angle=angle, frac=frac)
        if(missing(pos))
          pos <- if(abs(A[1] - B[1]) < abs(A[2] - B[2])) 4 else 3
      } else {
        ## two single-headed arrows with text 
        dM <- (shrink/2) * (B - A)
        AM <- M - dM
        BM <- M + dM
        newfrac <- frac/((1-shrink)/2)
        myarrows(AM[1], AM[2], A[1], A[2], ...,
               angle=angle, frac=newfrac, left=FALSE)
        myarrows(BM[1], BM[2], B[1], B[2], ...,
               angle=angle, frac=newfrac, left=FALSE)
      }
      text(M[1], M[2], attr(x, "txt"), cex=cex, pos=pos)
    }
    return(invisible(Window(x)))
  }
  plot.yardstick
})


print.yardstick <- function(x, ...) {
  splat("Yardstick")
  if(!is.null(txt <- attr(x, "txt")))
    splat("Text:", txt)
  ui <- summary(unitname(x))
  splat("Length:", pairdist(x)[1,2], ui$plural, ui$explain)
  splat("Midpoint:",
        paren(paste(signif(c(mean(x$x), mean(x$y)), 3), collapse=", ")))
  dx <- diff(range(x$x))
  dy <- diff(range(x$y))
  orient <- if(dx == 0) "vertical" else
            if(dy == 0) "horizontal" else
            paste(atan2(dy, dx) * 180/pi, "degrees")
  splat("Orientation:", orient)
  return(invisible(NULL))
}
