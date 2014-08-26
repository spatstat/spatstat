#
#	ppp.R
#
#	A class 'ppp' to define point patterns
#	observed in arbitrary windows in two dimensions.
#
#	$Revision: 4.98 $	$Date: 2014/07/29 02:59:13 $
#
#	A point pattern contains the following entries:	
#
#		$window:	an object of class 'owin'
#				defining the observation window
#
#		$n:	the number of points (for efficiency)
#	
#		$x:	
#		$y:	vectors of length n giving the Cartesian
#			coordinates of the points.
#
#	It may also contain the entry:	
#
#		$marks:	a vector of length n
#			whose entries are interpreted as the
#			'marks' attached to the corresponding points.	
#	
#--------------------------------------------------------------------------
ppp <- function(x, y, ..., window, marks, check=TRUE, drop=TRUE) {
  # Constructs an object of class 'ppp'
  #
  if(!missing(window))
    verifyclass(window, "owin")
  else
    window <- owin(...)

  if((missing(x) && missing(y)) || (length(x) == 0 && length(y) == 0))
    x <- y <- numeric(0)

  n <- length(x)
  if(length(y) != n)
    stop("coordinate vectors x and y are not of equal length")
  
  # validate x, y coordinates
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(y))
  ok <- is.finite(x) & is.finite(y)
  if(any(!ok)) {
    nbg <- is.na(x) | is.na(y)
    if(any(nbg)) {
      howmany <- if(all(nbg)) "all" else paste(sum(nbg),  "out of", length(nbg))
      stop(paste(howmany, "coordinate values are NA or NaN"))
    }
    howmany <- if(!any(ok)) "all" else paste(sum(!ok),  "out of", length(ok))
    stop(paste(howmany, "coordinate values are infinite"))
  }

  names(x) <- NULL
  names(y) <- NULL
  
  # check (x,y) points lie inside window
  if(check && n > 0) {
    ok <- inside.owin(x, y, window)
    nout <- sum(!ok)
    if(nout > 0) {
      warning(paste(nout,
                    ngettext(nout, "point was", "points were"),
                    "rejected as lying outside the specified window"))
      rr <- ripras(x,y)
      bb <- boundingbox(x,y)
      bb <- boundingbox(rr, bb, window)
      rejectwindow <-
        if(!is.null(rr)) rebound.owin(rr, bb) else bb
      rejects <- ppp(x[!ok], y[!ok], window=rejectwindow, check=FALSE)
      # discard illegal points
      x <- x[ok]
      y <- y[ok]
      n <- length(x)
    }
  } else nout <- 0
  # initialise ppp object
  pp <- list(window=window, n=n, x=x, y=y)
  # coerce marks to appropriate format
  if(missing(marks))
    marks <- NULL
  if(is.hyperframe(marks)) 
    stop("Hyperframes of marks are not implemented for ppp objects; use ppx")
  if(is.matrix(marks)) 
    marks <- as.data.frame(marks)
  ## drop dimensions?
  if(drop && is.data.frame(marks)) {
    nc <- ncol(marks)
    if(nc == 0)
      marks <- NULL
    else if(nc == 1)
      marks <- marks[,,drop=TRUE]
  }
  # attach marks 
  if(is.null(marks)) {
    # no marks
    pp$markformat <- "none"
  } else if(is.data.frame(marks)) {
    # data frame of marks
    pp$markformat <- "dataframe"
    if(nout > 0) {
      marks <- marks[ok, ]
      marks(rejects) <- marks[!ok,]
    }
    if(nrow(marks) != n)
      stop("number of rows of marks != length of x and y")
    pp$marks <- marks
  } else {
    # should be a vector or factor
    # To recognise vector, strip attributes
    isspecial <- is.factor(marks) ||
                 inherits(marks, "POSIXt") || inherits(marks, "Date")
    if(!isspecial)
      attributes(marks) <- NULL
    if(!(is.vector(marks) || isspecial))
      stop("Format of marks not understood")
    # OK, it's a vector or factor
    pp$markformat <- "vector"
    if(nout > 0) {
      marks(rejects) <- marks[!ok]
      marks <- marks[ok]
    }
    if(length(marks) != n)
      stop("length of marks vector != length of x and y")
    names(marks) <- NULL
    pp$marks <- marks
  }
  class(pp) <- "ppp"
  if(check && any(duplicated(pp)))
    warning("data contain duplicated points")
  if(nout > 0) 
    attr(pp, "rejects") <- rejects
  pp
}

#
#--------------------------------------------------------------------------
#

is.ppp <- function(x) { inherits(x, "ppp") }

#
#--------------------------------------------------------------------------
#

as.ppp <- function(X, ..., fatal=TRUE) {
  UseMethod("as.ppp")
}

as.ppp.ppp <- function(X, ..., fatal=TRUE) {
  check <- resolve.defaults(list(...), list(check=FALSE))$check
  return(ppp(X$x, X$y, window=X$window, marks=X$marks, check=check))
}

as.ppp.quad <- function(X, ..., fatal=TRUE) {
  return(union.quad(X))
}

as.ppp.data.frame <- function(X, W = NULL, ..., fatal=TRUE) {
  check <- resolve.defaults(list(...), list(check=TRUE))$check
  if(ncol(X) < 2) 
    return(complaining("X must have at least two columns",
                       fatal, value=NULL))

  if(is.null(W))
    return(complaining("x,y coords given but no window specified",
                       fatal, value=NULL))

  # columns 1 and 2 are assumed to be coordinates
  # marks from other columns
  marx <- if(ncol(X) > 2) X[, -(1:2)] else NULL

  if(is.function(W))
    Z <- cobble.xy(X[,1], X[,2], W, fatal, marks=marx, check=check)
  else {
    win <- as.owin(W)
    Z <- ppp(X[,1], X[,2], window = win, marks=marx, check=check)
  }

  return(Z)
}
    
as.ppp.matrix <- function(X, W = NULL, ..., fatal=TRUE) {
  check <- resolve.defaults(list(...), list(check=TRUE))$check
  if(!verifyclass(X, "matrix", fatal=fatal)
     || !is.numeric(X))
    return(complaining("X must be a numeric matrix",
                       fatal, value=NULL))

  if(ncol(X) < 2)
    return(complaining("X must have at least two columns",
                       fatal, value=NULL))

  if(is.null(W))
    return(complaining("x,y coords given but no window specified",
                       fatal, value=NULL))
    
  if(is.function(W))
    Z <- cobble.xy(X[,1], X[,2], W, fatal)
  else {
    win <- as.owin(W)
    Z <- ppp(X[,1], X[,2], window = win, check=check)
  }

  # add marks from other columns
  if(ncol(X) > 2)
    marks(Z) <- X[, -(1:2)]

  return(Z)
}
    
as.ppp.default <- function(X, W=NULL, ..., fatal=TRUE) {
	# tries to coerce data X to a point pattern
	# X may be:
	#	1. a structure with entries x, y, xl, xu, yl, yu
	#	2. a structure with entries x, y, area where
        #                    'area' has entries xl, xu, yl, yu
	#	3. a structure with entries x, y
        #       4. a vector of length 2, interpreted as a single point.
	# The second argument W is coerced to an object of class 'owin' by the 
	# function "as.owin" in window.S
        # If X also has an entry X$marks
        # then this will be interpreted as the marks vector for the pattern.
	#
  check <- resolve.defaults(list(...), list(check=TRUE))$check
  if(checkfields(X, c("x", "y", "xl", "xu", "yl", "yu"))) {
		xrange <- c(X$xl, X$xu)
		yrange <- c(X$yl, X$yu)
		if(is.null(X$marks))
			Z <- ppp(X$x, X$y, xrange, yrange, check=check)
		else
			Z <- ppp(X$x, X$y, xrange, yrange, 
				marks=X$marks, check=check)
		return(Z)
        } else if(checkfields(X, c("x", "y", "area"))
                  && checkfields(X$area, c("xl", "xu", "yl", "yu"))) {
                win <- as.owin(X$area)
                if (is.null(X$marks))
                  Z <- ppp(X$x, X$y, window=win, check=check)
                else
                  Z <- ppp(X$x, X$y, window=win, marks = X$marks, check=check)
                return(Z)
	} else if(checkfields(X, c("x", "y"))) {
                if(is.function(W))
                  return(cobble.xy(X$x, X$y, W, fatal))
		if(is.null(W)) {
                  if(fatal)
                    stop("x,y coords given but no window specified")
                  else
                    return(NULL)
                }
		win <- as.owin(W)
		if(is.null(X$marks))
                  Z <- ppp(X$x, X$y, window=win, check=check)
                else
                  Z <- ppp(X$x, X$y, window=win, marks=X$marks, check=check)
                return(Z)
        } else if(is.vector(X) && length(X) == 2) {
                win <- as.owin(W)
                Z <- ppp(X[1], X[2], window=win, check=check)
                return(Z)
	} else {
          if(fatal)
            stop("Can't interpret X as a point pattern")
          else
            return(NULL)
        }
}

cobble.xy <- function(x, y, f=ripras, fatal=TRUE, ...) {
  if(!is.function(f))
    stop("f is not a function")
  w <- f(x,y)
  if(!is.owin(w)) {
    gripe <- "Supplied function f did not return an owin object"
    if(fatal)
      stop(gripe)
    else {
      warning(gripe)
      return(NULL)
    }
  }
  return(ppp(x, y, window=w, ...))
}
  

# --------------------------------------------------------------

"[.ppp" <-
  function(x, i, j, drop, ...) {

        verifyclass(x, "ppp")
        
        if(missing(i) && missing(j))
          return(x)

        if(!missing(i)) {
          if(inherits(i, "owin")) {
            # i is a window
            window <- i
            ok <- inside.owin(x$x, x$y, window)
            x <- ppp(x$x[ok], x$y[ok], window=window, #SIC
                     marks=marksubset(x$marks, ok),
                     check=FALSE)
          } else if(inherits(i, "im")) {
            # i is an image
            if(i$type != "logical")
              stop(paste("Subset operator X[i] undefined",
                         "when i is a pixel image",
                         "unless it has logical values"), call.=FALSE)
            # convert logical image to window
            e <- sys.frame(sys.nframe())
            window <- solutionset(i, e)
            ok <- inside.owin(x$x, x$y, window)
            x <- ppp(x$x[ok], x$y[ok], window=window, #SIC
                     marks=marksubset(x$marks, ok),
                     check=FALSE)
          } else {
            # assume i is a subset index
            nx <- x$n
            if(nx == 0)
              return(x)
            subset <- seq_len(nx)[i]
            if(any(is.na(subset)))
              stop("Index out of bounds in [.ppp", call.=FALSE)
            x <- ppp(x$x[subset], x$y[subset], window=x$window,
                     marks=marksubset(x$marks, subset),
                     check=FALSE)
          } 
        }

        if(!missing(j))
          x <- x[j]   # invokes code above

        return(x)
}


# ------------------------------------------------------------------
#
#
scanpp <- function(filename, window, header=TRUE, dir="", multitype=FALSE) {
  filename <- if(dir=="") filename else
              paste(dir, filename, sep=.Platform$file.sep)
  df <- read.table(filename, header=header)
  if(header) {
    # check whether there are columns named 'x' and 'y'
    colnames <- dimnames(df)[[2]]
    xycolumns <- match(c("x", "y"), colnames, 0)
    named <- all(xycolumns > 0)
  } else {
    named <- FALSE
  }
  if(named) {
    x <- df$x
    y <- df$y
  } else {
    # assume x, y given in columns 1, 2 respectively
    x <- df[,1]
    y <- df[,2]
    xycolumns <- c(1,2)
  }
  if(ncol(df) == 2) 
      X <- ppp(x, y, window=window)
  else {
    marks <- df[ , -xycolumns]
    if(multitype) 
      marks <- factor(marks)
    X <- ppp(x, y, window=window, marks = marks)
  }
  X
}

#-------------------------------------------------------------------

"markspace.integral" <-
  function(X) {
  verifyclass(X, "ppp")
  if(!is.marked(X, dfok=TRUE))
    return(1)
  if(is.multitype(X))
    return(length(levels(marks(X))))
  else
    stop("Don't know how to compute total mass of mark space")
}

#-------------------------------------------------------------------

print.ppp <- local({

  graceful <- function(prefix, strings) {
    shortblurb <- paste(prefix, paste(strings, collapse=", "), "\n")
    if(nchar(shortblurb) < options("width")[[1]]) {
      cat(shortblurb)
    } else {
      cat(paste(prefix,"\n"))
      print(strings, quote=FALSE)
    }
    return(invisible(NULL))
  }
    
  print.ppp <- function(x, ...) {
    verifyclass(x, "ppp")
    ism <- is.marked(x, dfok=TRUE)
    cat(paste0(if(ism) "marked " else NULL,
               "planar point pattern:"),
        x$n,
        ngettext(x$n, "point", "points"),
        fill=TRUE)
    if(ism) {
      mks <- marks(x, dfok=TRUE)
      if(is.data.frame(mks)) {
        ## data frame of marks
        graceful("Mark variables:", names(mks))
      } else {
        ## vector of marks
        if(is.factor(mks)) {
          graceful("Multitype, with levels =", levels(mks))
        } else {
          ## Numeric, or could be dates
          if(inherits(mks, "Date")) {
            cat("marks are dates,",
                "of class", sQuote("Date"),
                fill=TRUE)
          } else if(inherits(mks, "POSIXt")) {
            cat("marks are dates,",
                "of class", sQuote("POSIXt"),
                fill=TRUE)
          } else {
            cat(paste0("marks are", if(is.numeric(mks)) " numeric," else NULL),
                "of storage type ", sQuote(typeof(mks)),
                fill=TRUE)
          }
        }
      }
    }
    print(x$window)
    if(!is.null(rejects <- attr(x, "rejects"))) {
      nrejects <- rejects$n
      cat("***",
          nrejects,
          ngettext(nrejects, "illegal point", "illegal points"),
          "stored in",
          paste("attr(,", dQuote("rejects"), ")", sep=""),
          "***",
          fill=TRUE)
    }
    if(!is.null(info <- attr(x, "info")) && inherits(info, "rmhInfoList"))
      cat("\nPattern was generated by",
          "Metropolis-Hastings simulation.",
          fill=TRUE)
    return(invisible(NULL))
  }

  print.ppp
})


summary.ppp <- function(object, ..., checkdup=TRUE) {
  verifyclass(object, "ppp")
  result <- list()
  result$is.marked <- is.marked(object, dfok=TRUE)
  result$n <- object$n
  result$window <- summary(object$window)
  result$intensity <- result$n/result$window$area
  if(checkdup) {
    result$nduplicated <- sum(duplicated(object))
    result$rounding <- rounding(object)
  }
  if(result$is.marked) {
    mks <- marks(object, dfok=TRUE)
    if(result$multiple.marks <- is.data.frame(mks)) {
      result$marknames <- names(mks)
      result$is.numeric <- FALSE
      result$marktype <- "dataframe"
      result$is.multitype <- FALSE
    } else {
      result$is.numeric <- is.numeric(mks)
      result$marknames <- "marks"
      result$marktype <- typeof(mks)
      result$is.multitype <- is.multitype(object)
    }
    if(result$is.multitype) {
      tm <- as.vector(table(mks))
      tfp <- data.frame(frequency=tm,
                        proportion=tm/sum(tm),
                        intensity=tm/result$window$area,
                        row.names=levels(mks))
      result$marks <- tfp
    } else 
      result$marks <- summary(mks)
  }
  class(result) <- "summary.ppp"
  if(!is.null(rejects <- attr(object, "rejects"))) 
    result$rejects <- rejects$n
  if(!is.null(info <- attr(object, "info")) && inherits(info, "rmhInfoList"))
    result$rmhinfo <- info
  return(result)
}

print.summary.ppp <- function(x, ..., dp=3) {
  verifyclass(x, "summary.ppp")
  terselevel <- spatstat.options("terse")
  splat(if(x$is.marked) "Marked planar" else "Planar",
        "point pattern: ",
        x$n,
        "points")
  oneline <- resolve.defaults(list(...), list(oneline=FALSE))$oneline
  if(oneline) return(invisible(NULL))
  unitinfo <- summary(x$window$units)
  splat("Average intensity",
        signif(x$intensity,dp),
        "points per square",
        unitinfo$singular,
        unitinfo$explain)
  ndup <- x$nduplicated
  if(waxlyrical('extras', terselevel) && !is.null(ndup) && (ndup > 0)) {
    parbreak(terselevel)
    splat("*Pattern contains duplicated points*")
  }
  rndg <- x$rounding
  if(waxlyrical('gory', terselevel) && !is.null(rndg)) {
    cat("\n")
    if(rndg >= 1) {
      cat("Coordinates are", "given to",
          rndg,
          "decimal", ngettext(rndg, "place", "places"),
          fill=TRUE)
      if(rndg <= 3) {
        cat("i.e. rounded to", "the nearest", "multiple of",
            10^(-rndg), unitinfo$plural, unitinfo$explain,
            fill=TRUE)
      }
    } else if(rndg == 0) {
      cat("Coordinates are", "integers", fill=TRUE)
      cat("i.e. rounded to", "the nearest", unitinfo$singular,
          unitinfo$explain, 
          fill=TRUE)
    } else {
      cat("Coordinates are", "multiples of",
          10^(-rndg), unitinfo$plural, unitinfo$explain, 
          fill=TRUE)
    }
    parbreak(terselevel)
  }
  if(x$is.marked) {
    if(x$multiple.marks) {
      splat("Mark variables:", commasep(x$marknames, ", "))
      cat("Summary:\n")
      print(x$marks)
    } else if(x$is.multitype) {
      cat("Multitype:\n")
      print(signif(x$marks,dp))
    } else {
      splat("marks are ",
            if(x$is.numeric) "numeric, ",
            "of type ", sQuote(x$marktype),
            sep="")
      cat("Summary:\n")
      print(x$marks)
    }
    parbreak(terselevel)
  }
  if(waxlyrical('extras', terselevel))
    print(x$window)
  if(waxlyrical('errors', terselevel) && !is.null(nrejects <- x$rejects)) {
    parbreak(terselevel)
    splat("***",
          nrejects,
          ngettext(nrejects, "illegal point", "illegal points"),
          "stored in",
          paste("attr(,", dQuote("rejects"), ")", sep=""),
          "***")
  }
  if(waxlyrical('gory', terselevel) && !is.null(info <- x$rmhinfo)) {
    cat("\nPattern was generated by",
        "Metropolis-Hastings algorithm rmh",
        fill=TRUE)
    print(info)
  }
  return(invisible(x))
}

# ---------------------------------------------------------------

identify.ppp <- function(x, ...) {
  verifyclass(x, "ppp")
  id <- identify(x$x, x$y, ...)
  if(!is.marked(x)) return(id)
  marks <- as.data.frame(x)[id, -(1:2)]
  out <- cbind(data.frame(id=id), marks)
  row.names(out) <- NULL
  return(out)
}

rebound <- function(x, rect) {
  UseMethod("rebound")
}

rebound.ppp <- function(x, rect) {
  verifyclass(x, "ppp")
  x$window <- rebound.owin(x$window, rect)
  return(x)
}

as.data.frame.ppp <- function(x, row.names=NULL, ...) {
  df <- data.frame(x=x$x, y=x$y, row.names=row.names)
  marx <- marks(x, dfok=TRUE)
  if(is.null(marx))
    return(df)
  if(is.data.frame(marx))
    df <- cbind(df, marx)
  else
    df <- data.frame(df, marks=marx)
  return(df)
}

is.empty.ppp <- function(x) { return(x$n == 0) }

npoints <- function(x) {
  UseMethod("npoints")
}

nobjects <- function(x) {
  UseMethod("nobjects")
}

nobjects.ppp <- npoints.ppp <- function(x) { x$n }


domain.ppp <- Window.ppp <- function(X, ...) { as.owin(X) }

"Window<-.ppp" <- function(X, ..., value) {
  verifyclass(value, "owin")
  return(X[value])
}

"Frame<-.ppp" <- function(X, value) {
  Frame(Window(X)) <- value
  return(X)
}

