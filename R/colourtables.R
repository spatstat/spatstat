#
# colourtables.R
#
# support for colour maps and other lookup tables
#
# $Revision: 1.26 $ $Date: 2013/07/14 07:08:20 $
#

colourmap <- function(col, ..., range=NULL, breaks=NULL, inputs=NULL) {
  # validate colour data 
  h <- col2hex(col)
  # store without conversion
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

plot.colourmap <- local({

  # recognised additional arguments to image.default() and axis()
  
  imageparams <- c("main", "asp", "sub", "axes", "ann",
                   "cex", "font", 
                   "cex.axis", "cex.lab", "cex.main", "cex.sub",
                   "col.axis", "col.lab", "col.main", "col.sub",
                   "font.axis", "font.lab", "font.main", "font.sub")
  axisparams <- c("cex", 
                  "cex.axis", "cex.lab",
                  "col.axis", "col.lab",
                  "font.axis", "font.lab")

  plot.colourmap <- function(x, ..., main,
                             xlim=NULL, ylim=NULL, vertical=FALSE, axis=TRUE,
                             labelmap=NULL) {
    if(missing(main))
      main <- short.deparse(substitute(x))
    stuff <- attr(x, "stuff")
    col <- stuff$outputs
    n   <- stuff$n
    discrete <- stuff$discrete
  #
    if(is.null(labelmap)) {
      labelmap <- function(x) x
    } else if(is.numeric(labelmap) && length(labelmap) == 1 && !discrete) {
      labscal <- labelmap
      labelmap <- function(x) { x * labscal }
    } else stopifnot(is.function(labelmap))
  
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
      # horizontal colour ribbon
      x <- linmap(v, rr, xlim)
      y <- ylim
      z <- matrix(v, ncol=1)
      do.call.matched("image.default",
                      resolve.defaults(list(x=x, y=y, z=z),
                                       list(...),
                                       list(main=main,
                                            xlim=xlim, ylim=ylim, asp=1.0,
                                            ylab="", xlab="", axes=FALSE,
                                            breaks=bks, col=col)),
                      extrargs=imageparams)
      if(axis) {
        # add horizontal axis
        if(discrete) {
          la <- paste(labelmap(stuff$inputs))
          at <- linmap(v, rr, xlim)
        } else {
          la <- prettyinside(rr)
          at <- linmap(la, rr, xlim)
          la <- labelmap(la)
        }
        # default axis position is below the ribbon (side=1)
        sidecode <- resolve.1.default("side", list(...), list(side=1))
        pos <- c(ylim[1], xlim[1], ylim[2], xlim[2])[sidecode]
        # draw axis
        do.call.matched("axis",
                        resolve.defaults(list(...),
                                         list(side = 1, pos = pos, at = at),
                                         list(labels=la)),
                        extrargs=axisparams)
      }
    } else {
      # vertical colour ribbon
      y <- linmap(v, rr, ylim)
      z <- matrix(v, nrow=1)
      x <- xlim
      do.call.matched("image.default",
                      resolve.defaults(list(x=x, y=y, z=z),
                                       list(...),
                                       list(main=main,
                                            ylim=ylim, xlim=xlim, asp=1.0,
                                            ylab="", xlab="", axes=FALSE,
                                            breaks=bks, col=col)),
                      extrargs=imageparams)
      if(axis) {
        # add vertical axis
        if(discrete) {
          la <- paste(labelmap(stuff$inputs))
          at <- linmap(v, rr, ylim)
        } else {
          la <- prettyinside(rr)
          at <- linmap(la, rr, ylim)
          la <- labelmap(la)
        }
        # default axis position is to the right of ribbon (side=4)
        sidecode <- resolve.1.default("side", list(...), list(side=4))
        pos <- c(ylim[1], xlim[1], ylim[2], xlim[2])[sidecode]
        # draw axis
        do.call.matched("axis",
                        resolve.defaults(list(...),
                                         list(side=4, pos=pos, at=at),
                                         list(labels=la)),
                        extrargs=axisparams)
      }
    }
    invisible(NULL)
  }

  plot.colourmap
})


# Interpolate a colourmap or lookup table defined on real numbers

interp.colourmap <- function(m, n=512) {
  if(!inherits(m, "colourmap"))
    stop("m should be a colourmap")
  st <- attr(m, "stuff")
  if(st$discrete) {
    # discrete set of input values mapped to colours
    xknots <- st$inputs
    # Ensure the inputs are real numbers
    if(!is.numeric(xknots))
      stop("Cannot interpolate: inputs are not numerical values")
  } else {
    # interval of real line, chopped into intervals, mapped to colours
    # Find midpoints of intervals
    bks <- st$breaks
    nb <- length(bks)
    xknots <- (bks[-1] + bks[-nb])/2
  }
  # corresponding colours in hsv coordinates
  yknots.hsv <- rgb2hsv(col2rgb(st$outputs))
  # transform 'hue' from polar to cartesian coordinate
  # divide domain into n equal intervals
  xrange <- range(xknots)
  xbreaks <- seq(xrange[1], xrange[2], length=n+1)
  xx <- (xbreaks[-1] + xbreaks[-(n+1)])/2
  # interpolate saturation and value in hsv coordinates
  yy.sat <- approx(x=xknots, y=yknots.hsv["s", ], xout=xx)$y
  yy.val <- approx(x=xknots, y=yknots.hsv["v", ], xout=xx)$y
  # interpolate hue by first transforming polar to cartesian coordinate
  yknots.hue <- 2 * pi * yknots.hsv["h", ]
  yy.huex <- approx(x=xknots, y=cos(yknots.hue), xout=xx)$y
  yy.huey <- approx(x=xknots, y=sin(yknots.hue), xout=xx)$y
  yy.hue <- (atan2(yy.huey, yy.huex)/(2 * pi)) %% 1
  # form colours using hue, sat, val
  yy <- hsv(yy.hue, yy.sat, yy.val)
  # done
  f <- colourmap(yy, breaks=xbreaks)
  return(f)
}

tweak.colourmap <- function(m, col, ..., inputs=NULL, range=NULL) {
  if(!inherits(m, "colourmap"))
    stop("m should be a colourmap")
  if(is.null(inputs) && is.null(range))
    stop("Specify either inputs or range")
  if(!is.null(inputs) && !is.null(range))
    stop("Do not specify both inputs and range")
  # determine indices of colours to be changed
  if(!is.null(inputs)) {
    ix <- m(inputs, what="index")
  } else {
    if(!(is.numeric(range) && length(range) == 2 && diff(range) > 0))
      stop("range should be a numeric vector of length 2 giving (min, max)")
    if(length(col2hex(col)) != 1)
      stop("When range is given, col should be a single colour value")
    ixr <- m(range, what="index")
    ix <- (ixr[1]):(ixr[2])
  }
  # reassign colours
  st <- attr(m, "stuff")
  outputs <- st$outputs
  is.hex <- function(z) identical(substr(z, 1, 7), substr(col2hex(z), 1, 7))
  result.hex <- FALSE
  if(is.hex(outputs)) {
    # convert replacement data to hex
    col <- col2hex(col)
    result.hex <- TRUE
  } else if(is.hex(col)) {
    # convert existing data to hex
    outputs <- col2hex(outputs)
    result.hex <- TRUE
  } else if(!(is.character(outputs) && is.character(col))) {
    # unrecognised format - convert both to hex
    outputs <- col2hex(outputs)
    col <- col2hex(col)
    result.hex <- TRUE
  }
  if(result.hex) {
    # hex codes may be 7 or 9 characters
    outlen <- nchar(outputs)
    collen <- nchar(col)
    if(length(unique(c(outlen, collen))) > 1) {
      # convert all to 9 characters
      if(any(bad <- (outlen == 7))) 
        outputs[bad] <- paste0(outputs[bad], "FF")
      if(any(bad <- (collen == 7))) 
        col[bad] <- paste0(col[bad], "FF")
    }
  }
  # Finally, replace
  outputs[ix] <- col
  st$outputs <- outputs
  attr(m, "stuff") <- st
  assign("stuff", st, envir=environment(m))
  return(m)
}

