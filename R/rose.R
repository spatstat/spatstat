#'
#'    rose.R
#'
#'   Rose diagrams
#'
#'   $Revision: 1.9 $  $Date: 2015/08/25 08:19:19 $
#'

rose <- function(x, ...) UseMethod("rose")

rose.default <- local({

  rose.default <- function(x, breaks = NULL, ...,
                           weights=NULL,
                           nclass=NULL,
                           unit=c("degree", "radian"),
                           start=0, clockwise=FALSE,
                           main) {
    if(missing(main) || is.null(main))
      main <- short.deparse(substitute(x))
    stopifnot(is.numeric(x))
    if(!is.null(weights))
      check.nvector(weights, length(x), things="observations")
    #' determine units
    missu <- missing(unit)
    unit <- match.arg(unit)
    unit <- validate.angles(x, unit, missu)
    FullCircle <- switch(unit, degree = 360, radian = 2*pi)
    #' reduce to [0, 2pi]
    x <- x %% FullCircle
    #' determine breakpoints strictly inside full circle
    breaks <- makebreaks(x, c(0, FullCircle), breaks, nclass)
    #' histogram without weights
    h <- do.call.matched(hist.default,
                         list(x=x, breaks=breaks, ..., plot=FALSE),
                         skipargs=graphicsAargh,
                         sieve=TRUE)
    result <- h$result
    otherargs <- h$otherargs
    #' redo weights, if given
    if(!is.null(weights)) {
      wh <- whist(x=x, breaks=breaks, weights=weights)
      result$count <- wh
      result$density <- wh/diff(breaks)
    }
    #
    do.call(rose.histogram,
            c(list(x=result, main=main,
                   unit=unit, start=start, clockwise=clockwise),
              otherargs))
  }

  graphicsAargh <- c("density", "angle", "col", "border",
                     "xlim", "ylim", "xlab", "ylab", "axes")

  makebreaks <- function(x, r, breaks=NULL, nclass=NULL) {
    use.br <- !is.null(breaks)
    if (use.br) {
      if (!is.null(nclass)) 
        warning("'nclass' not used when 'breaks' is specified")
    } else if (!is.null(nclass) && length(nclass) == 1L) {
      breaks <- nclass
    } else breaks <- "Sturges"
    use.br <- use.br && (nB <- length(breaks)) > 1L
    if (use.br) 
      breaks <- sort(breaks)
    else {
      if (is.character(breaks)) {
        breaks <- match.arg(tolower(breaks),
                            c("sturges", 
                              "fd",
                              "freedman-diaconis",
                              "scott"))
        breaks <- switch(breaks,
                         sturges = nclass.Sturges(x), 
                         `freedman-diaconis` = ,
                         fd = nclass.FD(x),
                         scott = nclass.scott(x), 
                         stop("unknown 'breaks' algorithm"))
      }
      else if (is.function(breaks)) {
        breaks <- breaks(x)
      }
      if (length(breaks) == 1) {
        if (!is.numeric(breaks) || !is.finite(breaks) || 
            breaks < 1L) 
          stop("invalid number of 'breaks'")
        breaks <- seq(r[1], r[2], length.out=breaks)
      }
      else {
        if (!is.numeric(breaks) || length(breaks) <= 1) 
          stop(gettextf("Invalid breakpoints produced by 'breaks(x)': %s", 
                        format(breaks)), domain = NA)
        breaks <- sort(breaks)
      }
    }
    return(breaks)
  }
  
  rose.default
})


rose.histogram <- function(x, ...,
                           unit=c("degree", "radian"),
                           start=0, clockwise=FALSE,
                           main, labels=TRUE, at=NULL, do.plot=TRUE) {
  if(missing(main) || is.null(main))
    main <- short.deparse(substitute(x))
  #' determine units
  missu <- missing(unit)
  unit <- match.arg(unit)
  #' validate
  bks <- x$breaks
  unit <- validate.angles(bks, unit, missu)
#  FullCircle <- switch(unit, degree = 360, radian = 2*pi)
  #' get sector sizes
  y <- x$density
  ymax <- max(y)
  #' draw disc
  insideclearance <- 0.1
  outsidespace <- if(!is.null(at) && length(at) == 0) 0 else
                  if(identical(labels, FALSE)) 0.1 else 0.25
  R <- (1+insideclearance) * ymax
  DD <- disc(R)
  Rout <- (1 + outsidespace) * R
  result <- do.call.matched(plot.owin,
                            resolve.defaults(list(x=disc(Rout),
                                                  main=main,
                                                  type="n"), 
                                             list(...)))
  do.call.matched(plot.owin,
                  resolve.defaults(list(x=DD,
                                        hatch=FALSE,
                                        add=TRUE),
                                   list(...)),
                  extrargs=graphicsPars("owin"),
                  skipargs="col")
  if(do.plot) {
    #' draw sectors
    ang <- ang2rad(bks, unit=unit, start=start, clockwise=clockwise)
    eps <- min(diff(ang), pi/128)/2
    for(i in seq_along(y)) {
      aa <- seq(ang[i], ang[i+1], by=eps)
      aa[length(aa)] <- ang[i+1]
      yi <- y[i]
      xx <- c(0, yi * cos(aa), 0)
      yy <- c(0, yi * sin(aa), 0)
      do.call.matched(polygon, list(x=xx, y=yy, ...))
    }
    #' add tick marks
    circticks(R, at=at, unit=unit, start=start, clockwise=clockwise,
              labels=labels)
  }
  #'
  return(invisible(result))
}

rose.density <- function(x, ..., unit=c("degree", "radian"),
                         start=0, clockwise=FALSE,
                         main, labels=TRUE, at=NULL, do.plot=TRUE) {
  if(missing(main) || is.null(main))
    main <- short.deparse(substitute(x))
  ang <- x$x
  rad <- x$y
  missu <- missing(unit)
  unit <- match.arg(unit)
  unit <- validate.angles(ang, unit, missu)
  #'
  result <- roseContinuous(ang, rad, unit, ...,
                           start=start, clockwise=clockwise,
                           main=main, labels=labels, at=at,
                           do.plot=do.plot)
  return(invisible(result))
}

rose.fv <- function(x, ..., unit=c("degree", "radian"),
                    start=0, clockwise=FALSE,
                    main, labels=TRUE, at=NULL, do.plot=TRUE) {
  if(missing(main) || is.null(main))
    main <- short.deparse(substitute(x))
  ang <- with(x, .x)
  rad <- with(x, .y)
  missu <- missing(unit)
  unit <- match.arg(unit)
  unit <- validate.angles(ang, unit, missu)
  #'
  result <- roseContinuous(ang, rad, unit, ...,
                           start=start, clockwise=clockwise,
                           main=main, labels=labels, at=at,
                           do.plot=do.plot)
  return(invisible(result))
}

roseContinuous <- function(ang, rad, unit, ...,
                           start=0, clockwise=FALSE,
                           main,
                           labels=TRUE, at=NULL,
                           do.plot=TRUE) {
  rmax <- max(rad)
  #' draw disc
  insideclearance <- 0.1
  outsidespace <- if(!is.null(at) && length(at) == 0) 0 else
                  if(identical(labels, FALSE)) 0.1 else 0.25
  R <- (1+insideclearance) * rmax
  DD <- disc(R)
  Rout <- (1 + outsidespace) * R
  result <- do.call.matched(plot.owin,
                            resolve.defaults(list(x=disc(Rout),
                                                  main=main,
                                                  type="n"), 
                                             list(...)))
  do.call.matched(plot.owin,
                  resolve.defaults(list(x=DD,
                                        add=TRUE,
                                        hatch=FALSE),
                                   list(...)),
                  extrargs=graphicsPars("owin"),
                  skipargs="col")
  #' draw plot
  if(do.plot) {
    ang <- ang2rad(ang, unit=unit, start=start, clockwise=clockwise)
    xx <- rad * cos(ang)
    yy <- rad * sin(ang)
    do.call.matched(polygon, list(x=xx, y=yy, ...), extrargs="lwd")
    circticks(R, at=at, unit=unit, start=start, clockwise=clockwise,
              labels=labels)
  }
  return(result)
}

ang2rad <- local({

  compasspoints <- c(E=0,N=90,W=180,S=270)
  
  ang2rad <- function(ang, unit=c("degree", "radian"),
                         start=0, clockwise=FALSE) {
    unit <- match.arg(unit)
    clocksign <- if(clockwise) -1 else 1
    stopifnot(length(start) == 1)
    if(is.character(start)) {
      if(is.na(match(toupper(start), names(compasspoints))))
        stop(paste("Unrecognised compass point", sQuote(start)), call.=FALSE)
      startdegrees <- compasspoints[[start]]
      start <- switch(unit,
                      degree = startdegrees,
                      radian = pi * (startdegrees/180))
      # start is measured anticlockwise
      ang <- start + clocksign * ang
    } else {
      stopifnot(is.numeric(start))
      # start is measured according to value of 'clockwise'
      ang <- clocksign * (start + ang)
    }
    rad <- switch(unit,
                  degree = pi * (ang/180),
                  radian = ang)
    return(rad)
  }

  ang2rad
})


circticks <- function(R, at=NULL, unit=c("degree", "radian"),
                      start=0, clockwise=FALSE, labels=TRUE) {
  unit <- match.arg(unit)
  FullCircle <- switch(unit, degree = 360, radian = 2*pi)
  if(is.null(at)) {
    at <- FullCircle * (0:23)/24
    major <- ((0:23) %% 6 == 0)
  } else {
    if(length(at) == 0) return(invisible(NULL))
    nat <- (at/FullCircle) * 4
    major <- abs(nat - round(nat)) < 0.01
  }
  atradians <- ang2rad(ang=at, unit=unit, start=start, clockwise=clockwise)
  tx <- R * cos(atradians)
  ty <- R * sin(atradians)
  expan <- ifelse(major, 1.1, 1.05)
  segments(tx, ty, expan * tx, expan * ty, lwd=major+1)
  if(!identical(labels, FALSE)) {
    if(identical(labels, TRUE)) {
      labels <- switch(unit,
                       degree=paste(round(at)),
                       radian=parse(text= simplenumber(at/pi, "pi", "*", 1e-3)))
    } else stopifnot(is.vector(labels) && length(labels) == length(at))
    big <- expan + 0.1
    text(big * tx, big * ty, labels=labels)
  }
  invisible(NULL)
}

validate.angles <- function(angles, unit=c("degree", "radian"), guess=TRUE) {
  #' validate
  width <- diff(range(angles))
  if(missing(unit) && guess && width <= 6.2832) {
    warning("Very small range of angles: treating them as radian")
    unit <- "radian"
  } else unit <- match.arg(unit)
  FullCircle <- switch(unit, degree = 360, radian = 2*pi)
  if(width > 1.002 * FullCircle)
    stop("Range of angles exceeds a full circle")
  return(unit)
}

