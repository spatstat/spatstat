#'
#'    rose.R
#'
#'   Rose diagrams
#'
#'   $Revision: 1.3 $  $Date: 2014/12/04 05:12:44 $
#'

rose <- function(x, ...) UseMethod("rose")

rose.default <- local({

  rose.default <- function(x, breaks = NULL, ...,
                           nclass=NULL,
                           unit=c("degree", "radian"), main) {
    if(missing(main) || is.null(main))
      main <- short.deparse(substitute(x))
    stopifnot(is.numeric(x))
    #' determine units
    missu <- missing(unit)
    unit <- match.arg(unit)
    unit <- validate.angles(x, unit, missu)
    FullCircle <- switch(unit, degree = 360, radian = 2*pi)
    #' reduce to [0, 2pi]
    x <- x %% FullCircle
    #' determine breakpoints strictly inside full circle
    breaks <- makebreaks(x, c(0, FullCircle), breaks, nclass)
    #'
    h <- do.call.matched(hist.default,
                         list(x=x, breaks=breaks, ..., plot=FALSE),
                         skipargs=graphicsAargh,
                         sieve=TRUE)
    do.call(rose.histogram, c(list(x=h$result, main=main), h$otherargs))
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
                           main, do.plot=TRUE) {
  if(missing(main) || is.null(main))
    main <- short.deparse(substitute(x))
  #' determine units
  missu <- missing(unit)
  unit <- match.arg(unit)
  #' validate
  bks <- x$breaks
  unit <- validate.angles(bks, unit, missu)
  FullCircle <- switch(unit, degree = 360, radian = 2*pi)
  #' get sector sizes
  y <- x$density
  ymax <- max(y)
  #' draw disc
  R <- 1.1 * ymax
  DD <- disc(R)
  result <- do.call.matched(plot.owin,
                            resolve.defaults(list(x=disc(R * 1.1),
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
    ang <- switch(unit, degree = pi * (bks/180), radian=bks)
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
    circticks(R)
  }
  #'
  return(invisible(result))
}

rose.density <- function(x, ..., unit=c("degree", "radian"),
                         main, do.plot=TRUE) {
  if(missing(main) || is.null(main))
    main <- short.deparse(substitute(x))
  ang <- x$x
  rad <- x$y
  missu <- missing(unit)
  unit <- match.arg(unit)
  unit <- validate.angles(ang, unit, missu)
  #'
  result <- roseContinuous(ang, rad, unit, ..., main=main, do.plot=do.plot)
  return(invisible(result))
}

rose.fv <- function(x, ..., unit=c("degree", "radian"),
                            main, do.plot=TRUE) {
  if(missing(main) || is.null(main))
    main <- short.deparse(substitute(x))
  ang <- with(x, .x)
  rad <- with(x, .y)
  missu <- missing(unit)
  unit <- match.arg(unit)
  unit <- validate.angles(ang, unit, missu)
  #'
  result <- roseContinuous(ang, rad, unit, ..., main=main, do.plot=do.plot)
  return(invisible(result))
}

roseContinuous <- function(ang, rad, unit, ..., main, do.plot=TRUE) {
  rmax <- max(rad)
  #' draw disc
  R <- 1.1 * rmax
  DD <- disc(R)
  result <- do.call.matched(plot.owin,
                            resolve.defaults(list(x=disc(R * 1.1),
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
    if(unit == "degree") ang <- pi * (ang/180)
    xx <- rad * cos(ang)
    yy <- rad * sin(ang)
    do.call.matched(polygon, list(x=xx, y=yy, ...), extrargs="lwd")
    circticks(R)
  }
  return(result)
}

circticks <- function(R, at) {
  if(missing(at)) {
    at <- 2 * pi * (0:23)/24
    major <- ((0:23) %% 6 == 0)
  } else {
    nat <- at * 2/pi
    major <- abs(nat - round(nat)) < 0.01
  }
  tx <- R * cos(at)
  ty <- R * sin(at)
  expan <- ifelse(major, 1.1, 1.05)
  segments(tx, ty, expan * tx, expan * ty, lwd=major+1)
  invisible(NULL)
}

validate.angles <- function(angles, unit, guess=TRUE) {
  #' validate
  width <- diff(range(angles))
  if(guess && width <= 6.2832) {
    warning("Very small range of angles: treating them as radian")
    unit <- "radian"
  }
  FullCircle <- switch(unit, degree = 360, radian = 2*pi)
  if(width > 1.002 * FullCircle)
    stop("Range of angles exceeds a full circle")
  return(unit)
}

