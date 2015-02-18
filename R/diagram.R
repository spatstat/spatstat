##
##   diagram.R
##
##   Simple objects for the elements of a diagram (text, arrows etc)
##    that are compatible with plot.layered and plot.listof
##
##   $Revision: 1.10 $ $Date: 2015/02/17 03:45:04 $

# ......... internal class 'diagramobj' supports other classes  .........

diagramobj <- function(X, ...) {
  if(inherits(try(Frame(X), silent=TRUE), "try-error"))
    stop("X is not a spatial object")
  a <- list(...)
  if(sum(nzchar(names(a))) != length(a))
    stop("All extra arguments must be named")
  attributes(X) <- append(attributes(X), a)
  class(X) <- c("diagramobj", class(X))
  return(X)
}

"[.diagramobj" <- function(x, ...) {
  y <- NextMethod("[")
  attributes(y) <- attributes(x)
  return(y)
}

shift.diagramobj <- function(X, ...) {
  y <- NextMethod("shift")
  attributes(y) <- attributes(X)
  return(y)
}

scalardilate.diagramobj <- function(X, f, ...) {
  y <- NextMethod("scalardilate")
  attributes(y) <- attributes(X)
  return(y)
}

# .............. user-accessible classes ................
# .........  (these only need a creator and a plot method) ......


## ...........  text .................

textstring <- function(x, y, txt=NULL, ...) {
  if(is.ppp(x) && missing(y)) {
    X <- x
    Window(X) <- boundingbox(x)
  } else {
    if(missing(y) && checkfields(x, c("x", "y"))) {
      y <- x$y
      x <- x$x
      stopifnot(length(x) == length(y))
    }
    X <- ppp(x, y, window=owin(range(x),range(y)))
  }
  marks(X) <- txt
  Y <- diagramobj(X, otherargs=list(...))
  class(Y) <- c("textstring", class(Y))
  return(Y)
}

plot.textstring <- function(x, ..., do.plot=TRUE) {
  txt <- marks(x)
  otha <- attr(x, "otherargs")
  if(do.plot) do.call.matched(text.default,
                              resolve.defaults(list(...),
                                               list(x=x$x, y=x$y, labels=txt),
                                               otha),
                              extrargs=c("srt", "family", "xpd"))
  return(invisible(Frame(x)))
}

print.textstring <- function(x, ...) {
  splat("Text string object")
  txt <- marks(x)
  if(npoints(x) == 1) {
    splat("Text:", dQuote(txt))
    splat("Coordinates:", paren(paste(as.vector(coords(x)), collapse=", ")))
  } else {
    splat("Text:")
    print(txt)
    splat("Coordinates:")
    print(coords(x))
  }
  return(invisible(NULL))
}
  
## ...........  'yardstick' to display scale information  ................

yardstick <- function(x0, y0, x1, y1, txt=NULL, ...) {
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
    X <- ppp(xx, yy, window=B, check=FALSE)
  }
  Window(X) <- boundingbox(X)
  Y <- diagramobj(X, txt=txt, otherargs=list(...))
  class(Y) <- c("yardstick", class(Y))
  return(Y)
}

plot.yardstick <- local({

  mysegments <- function(x0, y0, x1, y1, ..., moreargs=list()) {
    ## ignore unrecognised arguments without whingeing
    do.call.matched(segments,
                    resolve.defaults(list(x0=x0, y0=y0, x1=x1, y1=y1),
                                     list(...),
                                     moreargs),
                    extrargs=c("col", "lty", "lwd", "xpd", "lend"))
  }
  
  myarrows <- function(x0, y0, x1, y1, ...,
                       left=TRUE, right=TRUE,
                       angle=20, frac=0.25,
                       main, show.all, add) {
    mysegments(x0, y0, x1, y1, ...)
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
        mysegments(x0, y0, x0+up[1], y0+up[2], ...)
        mysegments(x0, y0, x0+lo[1], y0+lo[2], ...)
      }
      if(right) {
        mysegments(x1, y1, x1-up[1], y1-up[2], ...)
        mysegments(x1, y1, x1-lo[1], y1-lo[2], ...)
      }
    }
    return(invisible(NULL))
  }

  plot.yardstick <- function(x, ...,
                             angle=20,
                             frac=1/8,
                             split=FALSE,
                             shrink=1/4,
                             pos=NULL,
                             txt.args=list(),
                             txt.shift=c(0,0),
                             do.plot=TRUE) {
    if(do.plot) {
      txt <- attr(x, "txt")
      argh <- resolve.defaults(list(...), attr(x, "otherargs"))
      A <- as.numeric(coords(x)[1,])
      B <- as.numeric(coords(x)[2,])
      M <- (A+B)/2
      if(!split) {
        ## double-headed arrow
        myarrows(A[1], A[2], B[1], y1=B[2],
                 angle=angle, frac=frac, moreargs=argh)
        if(is.null(pos) && !("adj" %in% names(txt.args)))
          pos <- if(abs(A[1] - B[1]) < abs(A[2] - B[2])) 4 else 3
      } else {
        ## two single-headed arrows with text 
        dM <- (shrink/2) * (B - A)
        AM <- M - dM
        BM <- M + dM
        newfrac <- frac/((1-shrink)/2)
        myarrows(AM[1], AM[2], A[1], A[2],
                 angle=angle, frac=newfrac, left=FALSE, moreargs=argh)
        myarrows(BM[1], BM[2], B[1], B[2], 
                 angle=angle, frac=newfrac, left=FALSE, moreargs=argh)
      }
      if(is.null(txt.shift)) txt.shift <- rep(0, 2) else 
                             txt.shift <- ensure2vector(unlist(txt.shift))
      do.call.matched(text.default,
                      resolve.defaults(list(x=M[1] + txt.shift[1],
                                            y=M[2] + txt.shift[2]),
                                       txt.args,
                                       list(labels=txt, pos=pos),
                                       argh,
                                       .MatchNull=FALSE),
                      extrargs=c("srt", "family", "xpd"))
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


## code to draw a decent-looking arrow in spatstat diagrams
## (works in layered objects)

## The name 'onearrow' is used because R contains
## hidden functions [.arrow, length.arrow

onearrow <- function(x0, y0, x1, y1, txt=NULL, ...) {
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
    X <- ppp(xx, yy, window=B, check=FALSE)
  }
  Window(X) <- boundingbox(X)
  Y <- diagramobj(X, txt=txt, otherargs=list(...))
  class(Y) <- c("onearrow", class(Y))
  return(Y)
}

print.onearrow <- function(x, ...) {
  cat("Single arrow", fill=TRUE)
  if(!is.null(txt <- attr(x, "txt")))
    cat("Text:", txt, fill=TRUE)
  NextMethod("print")
}

plot.onearrow <- function(x, ...,
                          add=FALSE,
                          main="",
                          retract=0.05,   
                          headfraction=0.25,
                          headangle=12, # degrees
                          headnick=0.1, # fraction of head length
                          col.head=NA,
                          lwd.head=lwd,
                          lwd=1,
                          col=1,
                          zap=FALSE,
                          zapfraction=0.07,
                          pch=1, cex=1,
                          do.plot=TRUE,
                          do.points=FALSE,
                          show.all=!add) {
  result <- plot.ppp(x, main=main, add=add,
                     pch=pch, cex=cex,
                     do.plot=do.plot && do.points,
                     show.all=show.all)
  if(do.plot) {
    if(!do.points && !add)
      plot(Frame(x), main="", type="n")
    txt <- attr(x, "txt")
    argh <- resolve.defaults(list(...), attr(x, "otherargs"))
    A <- as.numeric(coords(x)[1,])
    B <- as.numeric(coords(x)[2,])
    V <- B - A
    AR <- A + retract * V
    BR <- B - retract * V
    H <- B - headfraction * V
    HN <- H + headnick * headfraction * V
    headlength <- headfraction * sqrt(sum(V^2))
    halfwidth <- headlength * tan((headangle/2) * pi/180)
    alpha <- atan2(V[2], V[1]) + pi/2
    U <- c(cos(alpha), sin(alpha))
    HL <- H + halfwidth * U
    HR <- H - halfwidth * U
    Head <- rbind(HN, HL, BR, HR, HN)
    if(!is.na(col.head))
      do.call.matched(polygon,
                      resolve.defaults(list(x=Head),
                                       argh,
                                       list(col=col.head, lwd=lwd.head)))
    if(!zap) {
      Tail <- AR
    } else {
      M <- (AR+HN)/2
      dM <- (zapfraction/2) * (1-headfraction) * V
      dM <- dM + c(-dM[2], dM[1])
      ML <- M + dM
      MR <- M - dM
      Tail <- rbind(AR, ML, MR)
    }
    do.call.matched(lines,
                    resolve.defaults(list(x=rbind(Tail, Head)),
                                     argh,
                                     list(col=col, lwd=lwd)),
                    extrargs=c("col", "lwd", "lty", "xpd", "lend"))
    if(!is.null(txt <- attr(x, "txt"))) {
      H <- (A+B)/2
      do.call.matched(text.default,
                      resolve.defaults(
                        list(x=H[1], y=H[2]),
                        argh,
                        list(labels=txt, pos=3 + (V[2] != 0))),
                      extrargs=c("srt", "family", "xpd"))
    }
  }
  return(invisible(result))
}
