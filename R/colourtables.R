#
# colourtables.R
#
# support for colour maps and other lookup tables
#
# $Revision: 1.16 $ $Date: 2011/05/18 01:34:03 $
#

colourmap <- function(col, ..., range=NULL, breaks=NULL, inputs=NULL) {
  # validate colour data 
  h <- col2hex(col)
  f <- lut(col, ..., range=range, breaks=breaks, inputs=inputs)
  class(f) <- c("colourmap", class(f))
  f
}

lut <- function(outputs, ..., range=NULL, breaks=NULL, inputs=NULL) {
  n <- length(outputs)
  given <- c(!is.null(range), !is.null(breaks), !is.null(inputs))
  names(given) <- c("range", "breaks", "inputs")
  ngiven <- sum(given)
  if(ngiven == 0)
    stop(paste("One of the arguments",
               sQuote("range"), ",", sQuote("breaks"), "or", sQuote("inputs"),
               "should be given"))
  if(ngiven > 1) {
    offending <- names(breaks)[given]
    stop(paste("The arguments",
               commasep(sQuote(offending)),
               "are incompatible"))
  }
  if(!is.null(inputs)) {
    # discrete set of input values mapped to output values
    stopifnot(length(inputs) == length(outputs))
    stuff <- list(n=n, discrete=TRUE, inputs=inputs, outputs=outputs)
    f <- function(x, what="value") {
      m <- match(x, stuff$inputs)
      if(what == "index")
        return(m)
      cout <- stuff$outputs[m]
      return(cout)
    }
  } else {
    # interval of real line mapped to colours
    if(is.null(breaks)) {
      breaks <- seq(from=range[1], to=range[2], length.out=length(outputs)+1)
    } else {
      stopifnot(is.numeric(breaks) && length(breaks) >= 2)
      stopifnot(length(breaks) == length(outputs) + 1)
      if(!all(diff(breaks) > 0))
        stop("breaks must be increasing")
    }
    stuff <- list(n=n, discrete=FALSE, breaks=breaks, outputs=outputs)
    f <- function(x, what="value") {
      stopifnot(is.numeric(x))
      x <- as.vector(x)
      z <- findInterval(x, stuff$breaks,
                        rightmost.closed=TRUE)
      if(what == "index")
        return(z)
      cout <- stuff$outputs[z]
      return(cout)
    }
  }
  attr(f, "stuff") <- stuff
  class(f) <- c("lut", class(f))
  f
}

print.lut <- function(x, ...) {
  stuff <- attr(x, "stuff")
  n <- stuff$n
  if(inherits(x, "colourmap")) {
    tablename <- "Colour map"
    outputname <- "colour"
  } else {
    tablename  <- "Lookup table"
    outputname <- "output"
  }
  if(stuff$discrete) {
    cat(paste(tablename, "for discrete set of input values\n"))
    out <- data.frame(input=stuff$inputs, output=stuff$outputs)
  } else {
    b <- stuff$breaks
    cat(paste(tablename, "for the range", prange(b[c(1,n+1)]), "\n"))
    leftend  <- rep("[", n)
    rightend <- c(rep(")", n-1), "]")
    inames <- paste(leftend, b[-(n+1)], ", ", b[-1], rightend, sep="")
    out <- data.frame(interval=inames, output=stuff$outputs)
  }
  colnames(out)[2] <- outputname
  print(out)
  invisible(NULL)
}

print.colourmap <- function(x, ...) {
  NextMethod("print")
}

summary.lut <- function(object, ...) {
  s <- attr(object, "stuff")
  if(inherits(object, "colourmap")) {
    s$tablename <- "Colour map"
    s$outputname <- "colour"
  } else {
    s$tablename  <- "Lookup table"
    s$outputname <- "output"
  }
  class(s) <- "summary.lut"
  return(s)
}

print.summary.lut <- function(x, ...) {
  n <- x$n
  if(x$discrete) {
    cat(paste(x$tablename, "for discrete set of input values\n"))
    out <- data.frame(input=x$inputs, output=x$outputs)
  } else {
    b <- x$breaks
    cat(paste(x$tablename, "for the range", prange(b[c(1,n+1)]), "\n"))
    leftend  <- rep("[", n)
    rightend <- c(rep(")", n-1), "]")
    inames <- paste(leftend, b[-(n+1)], ", ", b[-1], rightend, sep="")
    out <- data.frame(interval=inames, output=x$outputs)
  }
  colnames(out)[2] <- x$outputname
  print(out)  
}

plot.colourmap <- function(x, ..., main,
                           xlim=NULL, ylim=NULL, vertical=FALSE, axis=TRUE) {
  if(missing(main))
    main <- short.deparse(substitute(x))
  stuff <- attr(x, "stuff")
  col <- stuff$outputs
  n   <- stuff$n
  discrete <- stuff$discrete
  # determine pixel entries 'v' and colour map breakpoints 'bks'
  # to be passed to 'image.default'
  if(!discrete) {
    bks <- stuff$breaks
    rr <- range(bks)
    v <- seq(from=rr[1], to=rr[2], length.out=max(n+1, 1024))
  } else {
    v <- (1:n) - 0.5
    bks <- 0:n
    rr <- c(0,n)
  }
  # determine position of ribbon
  if(is.null(xlim) && is.null(ylim)) {
    u <- diff(rr)/10
    if(!vertical) {
      xlim <- rr
      ylim <- c(0,u)
    } else {
      xlim <- c(0,u)
      ylim <- rr
    }
  } else if(is.null(ylim)) {
    if(!vertical) 
      ylim <- c(0, diff(xlim)/10)
    else 
      ylim <- c(0, 10 * diff(xlim))
  } else if(is.null(xlim)) {
    if(!vertical) 
      xlim <- c(0, 10 * diff(ylim))
    else 
      xlim <- c(0, diff(ylim)/10)
  }
  # plot ribbon image
  linmap <- function(x, from, to) {
    to[1] + diff(to) * (x - from[1])/diff(from)
  }
  if(!vertical) {
    x <- linmap(v, rr, xlim)
    y <- ylim
    z <- matrix(v, ncol=1)
    do.call("image.default",
            resolve.defaults(list(x, y, z),
                             list(...),
                             list(main=main, xlim=xlim, ylim=ylim, asp=1.0,
                                  ylab="", xlab="", axes=FALSE,
                                  breaks=bks, col=col)))
    if(axis) {
      if(discrete) {
        la <- paste(stuff$inputs)
        at <- linmap(v, rr, xlim)
      } else {
        la <- pretty(rr)
        at <- linmap(la, rr, xlim)
      }
      axis(1, pos=ylim[1], at=at, labels=la)
    }
  } else {
    y <- linmap(v, rr, ylim)
    z <- matrix(v, nrow=1)
    x <- xlim
    do.call("image.default",
            resolve.defaults(list(x, y, z),
                             list(...),
                             list(main=main, ylim=ylim, xlim=xlim, asp=1.0,
                                  ylab="", xlab="", axes=FALSE,
                                  breaks=bks, col=col)))
    if(axis) {
      if(discrete) {
        la <- paste(stuff$inputs)
        at <- linmap(v, rr, ylim)
      } else {
        la <- pretty(rr)
        at <- linmap(la, rr, ylim)
      }
      axis(4, pos=xlim[2], at=at, labels=la)
    }
  }
  invisible(NULL)
}
