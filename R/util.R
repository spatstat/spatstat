#
#    util.S    miscellaneous utilities
#
#    $Revision: 1.187 $    $Date: 2015/10/01 10:23:50 $
#
#
matrowsum <- function(x) {
  x %*% rep.int(1, ncol(x))
}

matcolsum <- function(x) {
  rep.int(1, nrow(x)) %*% x
}
  
matrowany <- function(x) {
  (matrowsum(x) > 0)
}

matrowall <- function(x) {
  (matrowsum(x) == ncol(x))
}

matcolany <- function(x) {
  (matcolsum(x) > 0)
}

matcolall <- function(x) {
  (matcolsum(x) == nrow(x))
}

########
    # hm, this is SLOWER

apply23sum <- function(x) {
  dimx <- dim(x)
  if(length(dimx) != 3)
    stop("x is not a 3D array")
  result <- array(0, dimx[-1])

  nz <- dimx[3]
  for(k in 1:nz) {
    result[,k] <- matcolsum(x[,,k])
  }
  result
}
    
#######################
#
#    whist      weighted histogram
#

whist <- function(x, breaks, weights=NULL) {
    N <- length(breaks)
    if(length(x) == 0) 
      h <- numeric(N+1)
    else {
      # classify data into histogram cells (breaks need not span range of data)
      cell <- findInterval(x, breaks, rightmost.closed=TRUE)
      # values of 'cell' range from 0 to N.
      nb <- N + 1L
      if(is.null(weights)) {
        ## histogram
        h <- tabulate(cell+1L, nbins=nb)
      } else {
        ##  weighted histogram
        if(!spatstat.options("Cwhist")) {
          cell <- factor(cell, levels=0:N)
          h <- unlist(lapply(split(weights, cell), sum, na.rm=TRUE))
        } else {
          h <- .Call("Cwhist",
                     as.integer(cell), as.double(weights), as.integer(nb))
        }
      }
    }
    h <- as.numeric(h)
    y <- h[2:N]
    attr(y, "low") <- h[1]
    attr(y, "high") <- h[N+1]
    return(y)
}

######################
#
#   matrixsample         subsample or supersample a matrix
#

matrixsample <- function(mat, newdim, phase=c(0,0), scale, na.value=NA) {
  # 'phase+1' is the position of the [1,1] corner of the new matrix
  #  expressed in the coordinates of the old matrix.
  # 'scale' is the size of one step in the new matrix,
  #  expressed in the coordinates of the old matrix.
  # Both 'phase' and 'scale' can take any real value.
  olddim <- dim(mat)
  if(missing(scale)) scale <- (olddim - 1)/(newdim - 1)
  scale <- ensure2vector(scale)
  newdim  <- ensure2vector(newdim)
  newmat <- matrix(na.value, newdim[1], newdim[2])
  newrow <- 1:newdim[1]
  newcol <- 1:newdim[2]
  oldrow <- round(1 + phase[1] + (newrow-1) * scale[1])
  oldcol <- round(1 + phase[2] + (newcol-1) * scale[2])
  oldrow.ok <- (oldrow >= 1) & (oldrow <= olddim[1])
  oldcol.ok <- (oldcol >= 1) & (oldcol <= olddim[2])
  newmat[oldrow.ok, oldcol.ok] <- mat[oldrow[oldrow.ok],
                                      oldcol[oldcol.ok]]
  return(newmat)
}

# common invocation of matrixsample

rastersample <- function(X, Y) {
  stopifnot(is.im(X) || is.mask(X))
  stopifnot(is.im(Y) || is.mask(Y))
  phase <- c((Y$yrow[1] - X$yrow[1])/X$ystep,
             (Y$xcol[1] - X$xcol[1])/X$xstep)
  scale <- c(Y$ystep/X$ystep,
             Y$xstep/X$xstep)
  if(is.im(X)) {
    # resample an image
    if(!is.im(Y))
      Y <- as.im(Y)
    Xtype <- X$type
    Xv    <- X$v
    # handle factor-valued image as integer
    if(Xtype == "factor") 
      Xv <- array(as.integer(Xv), dim=X$dim)
    # resample
    naval <- switch(Xtype,
                 factor=,
                 integer= NA_integer_, 
                 logical = as.logical(NA_integer_), 
                 real = NA_real_, 
                 complex = NA_complex_, 
                 character = NA_character_,
                 NA)
    Y$v <- matrixsample(Xv, Y$dim, phase=phase, scale=scale, na.value=naval)
    # inherit pixel data type from X
    Y$type <- Xtype
    if(Xtype == "factor") {
      Y$v <- factor(Y$v, labels=levels(X))
      dim(Y$v) <- Y$dim
    }
  } else {
    # resample a mask
    if(!is.mask(Y)) Y <- as.mask(Y)
    Y$m <- matrixsample(X$m, Y$dim, phase=phase, scale=scale, na.value=FALSE)
  }
  return(Y)
}

pointgrid <- function(W, ngrid) {
  W <- as.owin(W)
  masque <- as.mask(W, dimyx=ngrid)
  rxy <- rasterxy.mask(masque, drop=TRUE)
  xx <- rxy$x
  yy <- rxy$y
  return(ppp(xx, yy, W))
}

# text magic

commasep <- function(x, join=" and ", flatten=TRUE) {
  px <- paste(x)
  nx <- length(px)
  if(nx <= 1) return(px)
  commas <- c(rep(", ", length(px)-2),
              join,
              "")
  out <- paste0(px, commas, collapse=if(flatten) "" else NULL)
  return(out)
}

paren <- function(x, type="(") {
  if(length(x) == 0) return(x)
  if(identical(type, "")) type <- "blank"
  switch(type,
         "(" = {
           out <- paste("(", x, ")", sep="")
         },
         "[" = {
           out <- paste("[", x, "]", sep="")
         },
         "{" = {
           out <- paste("{", x, "}", sep="")
         },
         blank = {
           out <- paste(x)
         },
         stop(paste("Unrecognised parenthesis type:", sQuote(type)))
         )
  out
}

unparen <- function(x) {
  x <- as.character(x)
  firstchar <- substr(x, 1, 1)
  n <- nchar(x)
  lastchar <- substr(x, n, n)
  enclosed <- n > 2 & (
                       (firstchar == "(" & lastchar == ")") |
                       (firstchar == "[" & lastchar == "]") |
                       (firstchar == "{" & lastchar == "}") )
  if(any(enclosed))
    x[enclosed] <- substr(x[enclosed], 2, n-1)
  return(x)
}

strsplitretain <- local({
  strsplitretain <- function(x, split=",") {
    ## split strings after occurrence of character b, but retain b
    y <- strsplit(x, split)
    lapply(y, addback, b=split)
  }
  addback <- function(x, b=",") {
    n <- length(x)
    if(n <= 1) x else c(paste0(x[-n], b), x[n])
  }    
  strsplitretain
})

truncline <- function(x, nc) {
  if(length(x) > 1)
    return(unlist(lapply(as.list(x), truncline, nc=nc)))
  ## split string into words
  y <- strsplit(x, " ", fixed=TRUE)[[1]]
  ## find max number of whole words that take up nc characters
  maxwords <- max(0, which(cumsum(nchar(y) + 1) <= nc+1))
  if(maxwords == length(y))
    return(x)
  ## truncation will occur.
  pad <- " [..]"
  nc <- nc - nchar(pad)
  maxwords <- max(0, which(cumsum(nchar(y) + 1) <= nc+1))
  z <- paste(y[seq_len(maxwords)], collapse=" ")
  d <- nc - nchar(z)
  if(d < 0)
    z <- substr(z, 1, nc)
  z <- paste0(z, pad)
  return(z)
}

padtowidth <- local({

  blankstring <- function(n) paste(rep(" ", n), collapse="")

  padtowidth <- function(a, b, justify=c("left", "right", "centre")) {
    justify <- match.arg(justify)
    if(is.character(b)) b <- nchar(b) else stopifnot(is.numeric(b))
    extra <- pmax(0, b - nchar(a))
    rpad <- lpad <- ""
    switch(justify,
           left = {
             rpad <- sapply(extra, blankstring)
           },
           right = {
             lpad <- sapply(extra, blankstring)
           },
           centre = {
             lpad <- sapply(floor(extra/2), blankstring)
             rpad <- sapply(ceiling(extra/2), blankstring)
           })
    result <- paste0(lpad, a, rpad)
    return(result)
  }

  padtowidth
})

fakecallstring <- function(fname, parlist) {
  cl <- do.call("call", append(list(name = fname), parlist))
  return(format(cl))
}

prange <- function(x) {
  stopifnot(length(x) == 2)
  paren(paste(x, collapse=", "), "[")
}

  
ordinal <- function(k) {
  last <- abs(k) %% 10
  lasttwo <- abs(k) %% 100
  isteen <- (lasttwo > 10 & lasttwo < 20)
  ending <- ifelse(isteen, "th",
                   ifelse(last == 1, "st",
                          ifelse(last == 2, "nd",
                                 ifelse(last == 3, "rd",
                                        "th"))))
  return(paste(k, ending, sep=""))
}

articlebeforenumber <- function(k) {
  k <- abs(k)
  if(k == 11) return("an")
  leading <- floor(k/10^floor(log10(k)))
  if(leading == 8) return("an")
  return("a")
}

# equivalent to rev(cumsum(rev(x)))

revcumsum <- function(x) {
  n <- length(x)
  if(identical(storage.mode(x), "integer")) {
    z <- .C("irevcumsum",
            x=as.integer(x),
            as.integer(n))
#            PACKAGE="spatstat")
    return(z$x)
  } else {
    z <- .C("drevcumsum",
            x=as.double(x),
            as.integer(n))
#            PACKAGE="spatstat")
    return(z$x)
  }
}

prolongseq <- function(x, newrange, step=NULL) {
  ## Extend a sequence x so that it covers the new range.
  stopifnot(length(newrange) == 2 && newrange[1] < newrange[2])
  ## Check 'x' is an evenly-spaced sequence
  if(length(x) > 1) {
    dx <- diff(x)
    if(any(dx <= 0))
      stop("x must be an increasing sequence")
    if(diff(range(dx)) > 0.01 * abs(mean(dx)))
      stop("x must be evenly spaced")
  }
  ## Infer step length
  if(!is.null(step)) {
    check.1.real(step)
    stopifnot(step > 0)
  } else if(length(x) > 1) {
    step <- mean(dx)
  } else stop("step is needed when x is a single value")

  ## 
  if(max(x) < newrange[1] || min(x) > newrange[2])
    stop("x lies entirely outside the desired range")
    
  ## add or trim data to left
  if(x[1] > newrange[1]) {
    leftbit <- seq(from=x[1], to=newrange[1], by= -step)[-1]
    x <- c(rev(leftbit), x)
    nleft <- length(leftbit)
  } else {
    nx <- length(x)
    x <- x[x >= newrange[1]]
    nleft <- length(x) - nx
  }

  # add or trim data to right
  nx <- length(x)
  if(newrange[2] > x[nx]) {
    rightbit <- seq(from=x[nx], to=newrange[2], by= step)[-1]
    x <- c(x, rightbit)
    nright <- length(rightbit)
  } else {
    x <- x[x <= newrange[2]]
    nright <- length(x) - nx
  }
  attr(x, "nleft") <- nleft
  attr(x, "nright") <- nright
  return(x)
}

## fill gaps in a sequence
fillseq <- function(x) {
  xname <- short.deparse(substitute(x))
  n <- length(x)
  if(n <= 1) return(x)
  rx <- range(x)
  dx <- diff(x)
  if(any(dx < 0)) stop(paste(xname, "should be an increasing sequence"),
                       call.=FALSE)
  ## guess step length
  eps <- diff(rx)/1e7
  step <- min(dx[dx > eps])
  ## make new sequence
  y <- seq(rx[1], rx[2], by=step)
  ny <- length(y)
  ## mapping from x to y
  i <- round((x - rx[1])/step) + 1L
  i <- pmin(ny, pmax(1, i))
  return(list(xnew=y, i=i))
}

intersect.ranges <- function(a, b, fatal=TRUE) {
  if(!is.null(a) && !is.null(b)) {
    lo <- max(a[1],b[1])
    hi <- min(a[2],b[2])
    if(lo <= hi)
      return(c(lo, hi))
  }
  if(fatal) stop("Intersection is empty")
  return(NULL)
}

inside.range <- function(x, r) {
  stopifnot(length(r) == 2 && r[1] <= r[2])
  return(x >= r[1] & x <= r[2])
}

check.in.range <- function(x, r, fatal=TRUE) {
  xname <- deparse(substitute(x))
  if(inside.range(x, r))
    return(TRUE)
  if(fatal) 
    stop(paste(xname, "should be a number between",
               r[1], "and", r[2]),
         call.=FALSE)
  return(FALSE)
}

startinrange <- function(x0, dx, r) {
  ## find y = x0 + n * dx such that y \in r
  if(all(inside.range(x0, r))) return(x0)
  stopifnot(is.numeric(dx) && length(dx) == 1)
  y <- x0 + dx * round((mean(r) - x0)/dx)
  y[!inside.range(y, r)] <- NA
  return(y)
}

prettyinside <- function(x, ...) {
  r <- range(x, na.rm=TRUE)
  if(diff(r) == 0) return(r[1])
  p <- pretty(x, ...)
  ok <- inside.range(p, r)
  return(p[ok])
}

prettydiscrete <- function(x, n=10) {
  nx <- length(x)
  dx <- nx %/% n
  if(dx < 1) return(x)
  i <- 1 + (0:(n-1)) * dx
  return(x[i])
}


check.range <- function(x, fatal=TRUE) {
  xname <- deparse(substitute(x))
  if(is.numeric(x) && identical(x, range(x, na.rm=TRUE)))
    return(TRUE)
  if(fatal) 
    stop(paste(xname, "should be a vector of length 2 giving (min, max)"))
  return(FALSE)
}
                        
niceround <- function(x, m=c(1,2,5,10)) {
  expo <- 10^as.integer(floor(log10(x)))
  y <- m * expo
  z <- y[which.min(abs(y - x))]
  return(z)
}

progressreport <- local({

  Put <- function(name, value, state) {
    if(is.null(state)) {
      putSpatstatVariable(paste0("Spatstat.", name), value)
    } else {
      state[[name]] <- value
    }
    return(state)
  }
  Get <- function(name, state) {
    if(is.null(state)) {
      value <- getSpatstatVariable(paste0("Spatstat.", name))
    } else {
      value <- state[[name]] 
    }
    return(value)
  }
  
  progressreport <- function(i, n, every=min(100,max(1, ceiling(n/100))),
                             nperline=min(charsperline,
                               every * ceiling(charsperline /(every+3))),
                             charsperline=60,
                             style=spatstat.options("progress"),
                             state=NULL) {
    missevery <- missing(every)
    if(i > n) {
      warning(paste("progressreport called with i =", i, "> n =", n))
      return(invisible(NULL))
    }
    if(style == "tk" && !requireNamespace("tcltk"))
      style <- "txtbar"
    switch(style,
           txtbar={
             if(i == 1) {
               ## initialise text bar
               state <- Put("ProgressBar",
                            txtProgressBar(1, n, 1, style=3),
                            state)
             } else {
               ## get text bar
               pbar <- Get("ProgressBar", state)
               ## update 
               setTxtProgressBar(pbar, i)
               if(i == n) {
                 close(pbar)
                 state <- Put("ProgressBar", NULL, state)
               } 
             }
           },
           tk={
             requireNamespace("tcltk")
             if(i == 1) {
               ## initialise text bar
               state <- Put("ProgressBar",
                            tcltk::tkProgressBar(title="progress",
                                                 min=0, max=n, width=300),
                            state)
             } else {
               ## get text bar
               pbar <- Get("ProgressBar", state)
               ## update 
               tcltk::setTkProgressBar(pbar, i,
                                       label=paste0(round(100 * i/n), "%"))
               if(i == n) {
                 close(pbar)
                 state <- Put("ProgressBar", NULL, state)
               } 
             }
           },
           tty={
             now <- proc.time()
             if(i == 1) {
               ## Initialise stuff
               if(missevery && every > 1 && n > 10) {
                 every <- niceround(every)
                 nperline <- min(charsperline,
                                 every * ceiling(charsperline /(every+3)))
               }
               showtime <- FALSE
               showevery <- n
               state <- Put("ProgressData",
                            list(every=every, nperline=nperline,
                                 starttime=now,
                                 showtime=FALSE, showevery=n),
                            state)
             } else {
               pd <- Get("ProgressData", state)
               if(is.null(pd))
                 stop(paste("progressreport called with i =", i,
                            "before i = 1"))
               every     <- pd$every
               nperline  <- pd$nperline
               showtime  <- pd$showtime
               showevery <- pd$showevery
               if(i < n) {
                 ## estimate time remaining
                 starttime <- pd$starttime
                 elapsed <- now - starttime
                 elapsed <- unname(elapsed[3])
                 rate <- elapsed/(i-1)
                 remaining <- rate * (n-i)
                 if(!showtime) {
                   ## show time remaining if..
                   if(rate > 20) {
                     ## .. rate is very slow
                     showtime <- TRUE
                     showevery <- 1
                   } else if(remaining > 180) {
                     ## ... more than 3 minutes remaining
                     showtime <- TRUE
                     showevery <- every
                     aminute <- ceiling(60/rate)
                     if(aminute < showevery) 
                       showevery <- min(niceround(aminute), showevery)
                   }
                 }
                 state <- Put("ProgressData",
                              list(every=every, nperline=nperline,
                                   starttime=starttime,
                                   showtime=showtime, showevery=showevery),
                              state)
               }
             }
             if(i == n) 
               cat(paste(" ", n, ".\n", sep=""))
             else if(every == 1 || i <= 3)
               cat(paste(i, ",", if(i %% nperline == 0) "\n" else " ", sep=""))
             else {
               if(i %% every == 0) 
                 cat(i)
               else
                 cat(".")
               if(i %% nperline == 0)
                 cat("\n")
             }
             if(i < n && showtime && (i %% showevery == 0)) {
               st <- paste("etd", codetime(round(remaining)))
               st <- paren(st, "[")
               cat(paste("", st, ""))
             }
             flush.console()
           },
           stop(paste("Unrecognised option for style:", dQuote(style)))
           )
    return(invisible(state))
  }

  progressreport
})
  
numalign <- function(i, nmax, zero="0") {
  stopifnot(i <= nmax)
  nplaces <- as.integer(ceiling(log10(nmax+1)))
  out <- paste(rep(zero, nplaces), collapse="")
  istring <- paste(i)
  ilen <- nchar(istring)
  substr(out, nplaces-ilen+1, nplaces) <- istring
  return(out)
}

as2vector <- function(x) {
  ## convert various wacky formats to numeric vector of length 2
  ## for use as coordinates of a single point.
  xname <- deparse(substitute(x))
  if(is.numeric(x)) {
    if(length(x) != 2)
      stop(paste(xname, "should have length 2"))
    return(x)
  }
  if(is.ppp(x)) {
    if(npoints(x) != 1)
      stop(paste(xname, "should consist of exactly one point"))
    return(c(x$x, x$y))
  }
  if(is.list(x) && all(c("x", "y") %in% names(x))) {
    if(length(x$x) != 1) stop(paste0(xname, "$x should have length 1"))
    if(length(x$y) != 1) stop(paste0(xname, "$y should have length 1"))
    return(c(x$x, x$y))
  }
  stop(paste("Format of", sQuote(xname), "not understood"))
}

ensure2vector <- function(x) {
  xname <- deparse(substitute(x))
  if(!is.numeric(x))
    stop(paste(xname, "is not numeric"))
  n <- length(x)
  if(n == 0 || n > 2)
    stop(paste(xname, "should be of length 1 or 2"))
  if(n == 1)
    return(rep(x,2))
  return(x)
}

ensure3Darray <- function(x) {
  nd <- length(dim(x))
  if(nd == 0) {
    x <- array(x, dim=c(length(x), 1, 1))
  } else if(nd == 2) {
    x <- array(x, dim=c(dim(x), 1))
  } else if(nd > 3) {
    laterdims <- dim(x)[-(1:3)]
    if(any(laterdims != 1))
      stop("Higher-dimensional array cannot be reduced to 3 dimensions")
    x <- array(x, dim=dim(x)[1:3])
  }
  return(x)
}

check.nvector <- function(v, npoints=NULL, fatal=TRUE, things="data points",
                          naok=FALSE, warn=FALSE, vname) {
  # vector of numeric values for each point/thing
  if(missing(vname))
    vname <- sQuote(deparse(substitute(v)))
  whinge <- NULL
  if(!is.numeric(v))
    whinge <- paste(vname, "is not numeric")
  else if(!is.atomic(v) || !is.null(dim(v)))  # vector with attributes
    whinge <- paste(vname, "is not a vector")
  else if(!is.null(npoints) && (length(v) != npoints))
    whinge <- paste("The length of", vname,
                    paren(paste0("=", length(v))), 
                    "should equal the number of", things,
                    paren(paste0("=", npoints)))
  else if(!naok && any(is.na(v)))
    whinge <- paste("Some values of", vname, "are NA or NaN")
  #
  if(!is.null(whinge)) {
    if(fatal) stop(whinge)
    if(warn) warning(whinge)
    ans <- FALSE
    attr(ans, "whinge") <- whinge
    return(ans)
  }
  return(TRUE)
}

check.nmatrix <- function(m, npoints=NULL, fatal=TRUE, things="data points",
                          naok=FALSE, squarematrix=TRUE, matchto="nrow",
                          warn=FALSE) {
  ## matrix of values for each thing or each pair of things
  mname <- sQuote(deparse(substitute(m)))
  whinge <- NULL
  if(!is.matrix(m))
    whinge <- paste(mname, "should be a matrix")
  else if(squarematrix && (nrow(m) != ncol(m)))
    whinge <- paste(mname, "should be a square matrix")
  else if(!naok && any(is.na(m)))
    whinge <- paste("Some values of", mname, "are NA or NaN")
  else if(!is.null(npoints)) {
    if(matchto=="nrow" && nrow(m) != npoints)
      whinge <- paste("Number of rows in", mname,
                      paren(paste0("=", nrow(m))),
                      "does not match number of", things,
                      paren(paste0("=", npoints)))
    else if(matchto=="ncol" && ncol(m) != npoints)
      whinge <- paste("Number of columns in", mname,
                      paren(paste0("=", ncol(m))),
                      "does not match number of", things,
                      paren(paste0("=", npoints)))
  }
  ##
  if(!is.null(whinge)) {
    if(fatal) stop(whinge)
    if(warn) warning(whinge)
    return(FALSE)
  }
  return(TRUE)
}

check.named.vector <- function(x, nam, context="", namopt=character(0)) {
  xtitle <- deparse(substitute(x))
  check.named.thing(x, nam, namopt, sQuote(xtitle),
                    is.numeric(x), "vector", context)
  opt <- namopt %in% names(x)
  return(x[c(nam, namopt[opt])])
}

check.named.list <- function(x, nam, context="", namopt=character(0)) {
  xtitle <- deparse(substitute(x))
  check.named.thing(x, nam, namopt, sQuote(xtitle),
                    is.list(x), "list", context)  
  opt <- namopt %in% names(x)
  return(x[c(nam, namopt[opt])])
}

check.named.thing <- function(x, nam, namopt=character(0), xtitle=NULL,
                              valid=TRUE, type="object", context="",
                              fatal=TRUE) {
  if(is.null(xtitle))
    xtitle <- sQuote(deparse(substitute(x)))
  # check whether names(x) contains all obligatory names 'nam'
  # and possibly some of the optional names 'namopt'
  namesx <- names(x)
  omitted <- !(nam %in% namesx)
  foreign <- !(namesx %in% c(nam, namopt))
  if(valid && !any(omitted) && !any(foreign))
    return(character(0))
  # some condition violated
  if(nzchar(context))
    xtitle <- paste(context, xtitle)
  whinge <- paste(xtitle,
                  "must be a named", paste(type, ",", sep=""),
                  "with components", commasep(nam))
  if(length(namopt) > 0)
    whinge <- paste(whinge, paren(paste("and optionally", commasep(namopt))))
  if(any(omitted)) {
    grizzle <- paste(ngettext(sum(omitted), "parameter", "parameters"),
                     commasep(nam[omitted]),
                     "omitted")
    whinge <- paste(whinge, grizzle, sep="; ")
  }
  if(any(foreign)) {
    grizzle <- paste(ngettext(sum(foreign), "component", "components"),
                     commasep(namesx[foreign]),
                     "not recognised")
    whinge <- paste(whinge, grizzle, sep="; ")
  }
  if(fatal)
    stop(whinge, call.=FALSE)
  return(whinge)
}


forbidNA <- function(x, context="", xname, fatal=TRUE, usergiven=TRUE) {
  if(missing(xname)) xname <- sQuote(deparse(substitute(x)))
  if(any(is.na(x))) {
    if(usergiven) {
      # argument came from user
      offence <- ngettext(length(x), "be NA", "contain NA values")
      whinge <- paste(context, xname, "must not", offence)
    } else {
      # argument was computed internally
      violates <- ngettext(length(x), "is NA", "contains NA values")
      whinge <- paste(context, xname, violates)
    }
    if(fatal) stop(whinge, call.=FALSE)
    warning(whinge, call.=FALSE)
    return(FALSE)
  }
  return(TRUE)
}

check.finite <- function(x, context="", xname, fatal=TRUE, usergiven=TRUE) {
  if(missing(xname)) xname <- sQuote(deparse(substitute(x)))
  forbidNA(x, context, xname, fatal=fatal, usergiven=usergiven)
  if(any(!is.finite(x))) {
    if(usergiven) {
      # argument came from user
      oblige <- ngettext(length(x),
                         "be a finite value", "contain finite values")
      whinge <- paste(context, xname, "must", oblige)
    } else {
      # argument was computed internally
      violates <- ngettext(length(x),
                           "is not finite", "contains non-finite values")
      whinge <- paste(context, xname, violates)
    }
    if(fatal) stop(whinge, call.=FALSE)
    warning(whinge, call.=FALSE)
    return(FALSE)
  }
  return(TRUE)
}

evenly.spaced <- function(x, tol=1e-07) {
  # test whether x is evenly spaced and increasing
  dx <- diff(x)
  if(any(dx <= .Machine$double.eps))
    return(FALSE)
  # The following test for equal spacing is used in hist.default
  if(diff(range(dx)) > tol * mean(dx))
    return(FALSE)
  return(TRUE)
}

adjustthinrange <- function(ur,vstep,vr) {
  if(diff(ur) >= vstep) return(ur)
  ur <- mean(ur) + c(-1,1) * vstep/2
  if(ur[1] < vr[1]) ur <- vr[1] + c(0,1)*vstep
  if(ur[2] > vr[2]) ur <- vr[2] - c(1,0)*vstep
  return(ur)
}

validposint <- function(n, caller, fatal=TRUE) {
  nname <- deparse(substitute(n))
  if(length(n) != 1 || n != round(n) || n <=0) {
    if(!fatal)
      return(FALSE)
    prefix <- if(!missing(caller)) paste("In ", caller, ",", sep="") else NULL
    stop(paste(prefix, nname, "should be a single positive integer"),
         call.=FALSE)
  }
  return(TRUE)
}

# wrangle data.frames

findfirstfactor <- function(x) {
  stopifnot(is.data.frame(x) || is.hyperframe(x))
  isfac <- unlist(lapply(as.list(x), is.factor))
  if(!any(isfac)) 
    return(NULL)
  min(which(isfac))
}

firstfactor <- function(x) {
  j <- findfirstfactor(x)
  if(is.null(j)) return(NULL)
  return(x[, j, drop=TRUE])
}

onecolumn <- function(m) {
  switch(markformat(m),
         none=stop("No marks provided"),
         vector=m,
         dataframe=m[,1, drop=TRUE],
         NA)
}

# errors and checks

complaining <- function(whinge, fatal=FALSE, value=NULL) {
  if(fatal) stop(whinge, call.=FALSE)
  warning(whinge, call.=FALSE)
  return(value)
}

check.1.real <- function(x, context="", fatal=TRUE) {
  xname <- deparse(substitute(x))
  if(!is.numeric(x) || length(x) != 1) {
    whinge <-  paste(sQuote(xname), "should be a single number")
    if(nzchar(context)) whinge <- paste(context, whinge)
    return(complaining(whinge, fatal=fatal, value=FALSE))
  }
  return(TRUE)
}

check.1.integer <- function(x, context="", fatal=TRUE) {
  xname <- deparse(substitute(x))
  if(!is.numeric(x) || length(x) != 1 || !is.finite(x) || x %% 1 != 0) {
    whinge <-  paste(sQuote(xname), "should be a single finite integer")
    if(nzchar(context)) whinge <- paste(context, whinge)
    return(complaining(whinge, fatal=fatal, value=FALSE))
  }
  return(TRUE)
}

explain.ifnot <- function(expr, context="") {
  ex <- deparse(substitute(expr))
  ans <- expr
  if(!(is.logical(ans) && length(ans) == 1 && ans))
    stop(paste(context, "it must be TRUE that", sQuote(ex)), call.=FALSE)
}

warn.ignored.args <- function(..., context=NULL) {
  if((narg <- length(list(...))) > 0) {
    whinge <- paste(narg, "unrecognised",
                    ngettext(narg, "argument was", "arguments were"),
                    "ignored")
    if(!is.null(context)) whinge <- paste(context, whinge)
    warning(context)
  }
}

multiply.only.finite.entries <- function(x, a) {
  # In ppm a potential value that is -Inf must remain -Inf
  # and a potential value that is 0 multiplied by NA remains 0
  y <- x
  ok <- is.finite(x) & (x != 0)
  y[ok] <- a * x[ok]
  return(y)
}

singlestring <- function(s, coll="") {
  s <- as.character(s)
  if(length(s) > 1)
    s <- paste(s, collapse=coll)
  return(s)
}

verbalogic <- function(x, op="and") {
  stopifnot(is.character(x))
  istrue <- (x == "TRUE")
  isfalse <- (x == "FALSE")
  isvariable <- !istrue & !isfalse
  y <- x[isvariable]
  switch(op,
         and={
           if(any(isfalse))
             return("FALSE")
           if(all(istrue))
             return("TRUE")
           return(paste(y, collapse=" and "))
         },
         or={
           if(all(isfalse))
             return("FALSE")
           if(any(istrue))
             return("TRUE")
           return(paste(y, collapse=" or "))
         },
         not={
           x[isfalse] <- "TRUE"
           x[istrue] <- "FALSE"
           x[isvariable] <- paste("not {", y, "}")
         },
         stop(paste("Unrecognised operation", sQuote(op))))
}

sensiblevarname <- function(guess, fallback, maxlen=12) {
  out <- if(is.character(guess) &&
            length(guess) == 1  &&
            make.names(guess) == guess) guess else fallback
  out <- substr(out, 1, maxlen)
  return(out)
}

## deparse() can sometimes be equivalent to dumping the whole object
short.deparse <- function(x, maxlen=60) {
  deparse(x,
          nlines=1,
          width.cutoff=maxlen,
          control="delayPromises")
}

## deparse() can produce multiple lines of text
flat.deparse <- function(x) {
  y <- paste(deparse(x), collapse=" ")
  y <- gsub("\n", " ", y)
  y <- gsub(" ", "", y)
  return(y)
}

good.names <- function(nama, defaults, suffices) {
  # ensure sensible, unique names 
  stopifnot(is.character(defaults))
  if(!missing(suffices))
    defaults <- paste(defaults, suffices, sep="")
  result <- nama
  if(is.null(result))
    result <- defaults
  else if(any(blank <- !nzchar(result)))
    result[blank] <- defaults[blank]
  if(anyDuplicated(result))
    result <- make.names(result, unique=TRUE)
  return(result)
}

cat.factor <- function (..., recursive=FALSE) {
  lll <- list(...)
  chk <- sapply(lll,is.factor)
  if(!all(chk))
    stop("First argument is a factor and at least one other argument is not.\n")
  lll <- lapply(lll,as.data.frame,nm="v1")
  return(do.call(rbind,lll)[,1])
}

nzpaste <- function(..., sep=" ", collapse=NULL) {
  # Paste only the non-empty strings
  v <- list(...)
  ok <- unlist(lapply(v, function(z) {any(nzchar(z))}))
  do.call("paste", append(v[ok], list(sep=sep, collapse=collapse)))
}

substringcount <- function(x, y) {
  ## count occurrences of 'x' in 'y'
  yy <- paste0("a", y, "a")
  splot <- strsplit(yy, split=x, fixed=TRUE)
  nhits <- unlist(lapply(splot, length)) - 1
  return(nhits)
}

is.parseable <- function(x) {
  unlist(lapply(x, function(z) {
    !inherits(try(parse(text=z), silent=TRUE), "try-error")
  }))
}

make.parseable <- function(x) {
  if(all(is.parseable(x))) x else make.names(x)
}

# paste(expression(..)) seems to be broken

paste.expr <- function(x) {
  unlist(lapply(x, function(z) { paste(deparse(z), collapse="") }))
}

#   gsub(".", replacement, x) but only when "." appears as a variable

gsubdot <- function(replacement, x) {
  x <- as.character(x)
  stopifnot(length(x) == 1)
  # find all positions of "." in x
  dotpos <- gregexpr("\\.", x)[[1]]
  if(all(dotpos == -1)) return(x)
  # find all positions of "." preceded or followed by alphanumeric
  dotbefore <- gregexpr("\\.[0-9A-Za-z]", x)[[1]]
  dotafter <- gregexpr("[0-9A-Za-z]\\.", x)[[1]] - 1
  # exclude them
  dotpos <- setdiff(dotpos, union(dotbefore, dotafter))
  #
  if(length(dotpos) == 0) return(x)
  lenrep <-length(replacement)
  while(length(dotpos) > 0) {
    dp <- dotpos[1]
    x <- paste0(substr(x, 0, dp-1), replacement, substr(x, dp+1, nchar(x)))
    dotpos <- dotpos[-1] + lenrep-1
  }
  return(x)
}

badprobability <- function(x, NAvalue=NA) {
  ifelse(is.na(x), NAvalue, !is.finite(x) | x < 0 | x > 1)
}

# test for equivalence of two functions 
samefunction <- function(f, g) {
  identical(deparse(f), deparse(g))
}

codetime <- local({
  uname <- c("min", "hours", "days", "years",
             "thousand years", "million years", "billion years")
  u1name <- c("min", "hour", "day", "year",
             "thousand years", "million years", "billion years")
  multiple <- c(60, 60, 24, 365, 1e3, 1e3, 1e3)
  codehms <- function(x) {
    sgn <- if(x < 0) "-" else ""
    x <- round(abs(x))
    hours <- x %/% 3600
    mins  <- (x %/% 60) %% 60
    secs  <- x %% 60
    h <- if(hours > 0) paste(hours, ":", sep="") else ""
    started <- (hours > 0)
    m <- if(mins > 0) {
      paste(if(mins < 10 && started) "0" else "", mins, ":", sep="")
    } else if(started) "00:" else ""
    started <- started | (mins > 0)
    s <- if(secs > 0) {
      paste(if(secs < 10 && started) "0" else "", secs, sep="")
    } else if(started) "00" else "0"
    if(!started) s <- paste(s, "sec")
    paste(sgn, h, m, s, sep="")
  }
  codetime <- function(x, hms=TRUE, what=c("elapsed","user","system")) {
    if(inherits(x, "proc_time")) x <- summary(x)[[match.arg(what)]] 
    if(!is.numeric(x) || length(x) != 1)
      stop("codetime: x must be a proc_time object or a single number")
    sgn <- if(x < 0) "-" else ""
    x <- abs(x)
    if(x < 60)
      return(paste(sgn, signif(x, 3), " sec", sep=""))
    # more than 1 minute: round to whole number of seconds
    x <- round(x)
    if(hms && (x < 60 * 60 * 24))
      return(paste(sgn, codehms(x), sep=""))
    u <- u1 <- "sec"
    for(k in seq_along(multiple)) {
      if(x >= multiple[k]) {
        x <- x/multiple[k]
        u <- uname[k]
        u1 <- u1name[k]
      } else break
    }
    xx <- round(x, 1)
    ux <- if(xx == 1) u1 else u
    paste(sgn, xx, " ", ux, sep="")
  }
  codetime
})

# defines the current favorite algorithm for 'order' 
fave.order <- function(x) { sort.list(x, method="quick", na.last=NA) }

# convert any appropriate subset index for a point pattern
# to a logical vector

ppsubset <- function(X, I) {
  Iname <- deparse(substitute(I))
  # I could be a function to be applied to X
  if(is.function(I)) {
    I <- I(X)
    if(!is.vector(I)) {
      warning(paste("Function", sQuote(Iname), "did not return a vector"),
              call.=FALSE)
      return(NULL)
    }
  }      
  # I is now an index vector
  n <- npoints(X)
  i <- try(seq_len(n)[I])
  if(inherits(i, "try-error") || any(is.na(i))) {
    warning(paste("Invalid subset index", sQuote(Iname)),
            call.=FALSE)
    return(NULL)
  }
  if(is.logical(I))
    return(I)
  # convert to logical
  Z <- rep.int(FALSE, n)
  Z[I] <- TRUE
  return(Z)
}


trap.extra.arguments <- function(..., .Context="", .Fatal=FALSE) {
  z <- list(...)
  if((narg <- length(z)) == 0) return(FALSE)
  nama <- names(z)
  named <- nzchar(nama)
  whinge <- paste(.Context, ":", sep="")
  if(any(named)) {
    # some arguments are named: ensure all are named
    nama <- sQuote(nama)
    if(!all(named)) 
      nama[!named] <- paste("[Arg", 1:length(nama), ,"]", sep="")[!named]
    whinge <- paste(whinge,
                    "unrecognised",
                    ngettext(narg, "argument", "arguments"),
                    commasep(nama),
                    ngettext(narg, "was", "were"), "ignored")
  } else {
    # all arguments unnamed
    whinge <- paste(whinge, 
                    narg, "unrecognised",
                    ngettext(narg, "argument was", "arguments were"),
                    "ignored")   
  }
  if(.Fatal) stop(whinge, call.=FALSE) else warning(whinge, call.=FALSE)
  return(TRUE)
}

dotexpr.to.call <- function(expr, dot="funX", evaluator="eval.fv") {
  # convert an expression into a function call
  # replacing "." by the specified variable 
  stopifnot(is.expression(expr))
  aa <- substitute(substitute(ee, list(.=as.name(d))),
                   list(ee=expr, d=dot))
  bb <- eval(parse(text=deparse(aa)))
  cc <- as.call(bb)
  cc[[1]] <- as.name("eval.fv")
  return(cc)
}

inject.expr <- function(base, expr) {
  ## insert an expression inside a call and parse it
  txt <- sub(".", as.character(expr), as.character(base), fixed=TRUE)
  parse(text=txt)
}

  
## Match variable names to objects in 'data' list or environment
getdataobjects <- function(nama, envir, datalist=NULL, fatal=FALSE) {
  if(is.null(nama)) return(NULL)
  stopifnot(is.character(nama))
  n <- length(nama)
  y <- vector(mode="list", length=n)
  names(y) <- nama
  if(!is.null(datalist)) {
    hit <- nama %in% names(datalist)
    if(any(hit))
      y[hit] <- as.list(datalist)[nama[hit]]
    external <- unlist(lapply(y, is.null))
  } else external <- rep(TRUE, n)
  y[external] <- mget(nama[external], envir=envir,
                    ifnotfound=list(NULL), inherits=TRUE)
  if(fatal && any(bad <- unlist(lapply(y, is.null)))) {
    nbad <- sum(bad)
    stop(paste(ngettext(nbad, "Covariate", "Covariates"),
               commasep(sQuote(nama[bad])),
               ngettext(nbad, "was not found", "were not found")),
         call.=FALSE)
  }
  names(y) <- nama
  attr(y, "external") <- external
  return(y)
}
 
## print names and version numbers of libraries loaded

sessionLibs <- function() {
  a <- sessionInfo()
  b <- unlist(lapply(a$otherPkgs, getElement, name="Version"))
  g <- rbind(names(b), unname(b))
  d <- apply(g, 2, paste, collapse=" ")
  if(length(d) > 0) {
    cat("Libraries loaded:\n")
    for(di in d) cat(paste("\t", di, "\n"))
  } else cat("Libraries loaded: None\n")
  return(invisible(d))
}

dropifsingle <- function(x) if(length(x) == 1) x[[1]] else x

# timed objects

timed <- function(x, ..., starttime=NULL, timetaken=NULL) {
  if(is.null(starttime)) # time starts now.
    starttime <- proc.time()
  # evaluate expression if any
  object <- x
  timetaken <- proc.time() - starttime
  class(object) <- c("timed", class(object))
  attr(object, "timetaken") <- timetaken
  return(object)
}

print.timed <- function(x, ...) {
  # strip the timing information and print the rest.
  taken <- attr(x, "timetaken")
  cx <- class(x)
  attr(x, "timetaken") <- NULL
  class(x) <- cx[cx != "timed"]
  NextMethod("print")
  # Now print the timing info
  cat(paste("\nTime taken:", codetime(taken), "\n"))
  return(invisible(NULL))
}

# wrapper for computing weighted variance of a vector
# Note: this includes a factor 1 - sum(v^2) in the denominator
# where v = w/sum(w). See help(cov.wt)

weighted.var <- function(x, w, na.rm=FALSE) {
  bad <- is.na(w) | is.na(x)
  if(any(bad)) {
    if(!na.rm) return(NA_real_)
    ok <- !bad
    x <- x[ok]
    w <- w[ok]
  }
  cov.wt(matrix(x, ncol=1),w)$cov[]
}

# efficient replacements for ifelse()
# 'a' and 'b' are single values
# 'x' and 'y' are vectors of the same length as 'test'

# ifelse(test, a, b)
ifelseAB <- function(test,  a, b) {
  y <- rep.int(b, length(test))
  y[test] <- a
  return(y)
}

# ifelse(test, a, x)
ifelseAX <- function(test, a, x) {
  y <- x
  y[test] <- a
  return(y)
}

# ifelse(test, x, b)
ifelseXB <- function(test, x, b) {
  y <- rep.int(b, length(test))
  y[test] <- x[test]
  return(y)
}
  
# ifelse(test, x, y)
ifelseXY <- function(test, x, y) {
  z <- y
  z[test] <- x[test]
  return(z)
}

#.... very special cases ......

# ifelse(test, 1, NA)
ifelse1NA <- function(test) {
  y <- as.integer(test)
  y[!test] <- NA
  return(y)
}

# ifelse(test, 0, NA)
ifelse0NA <- function(test) {
  nyet <- !test
  y <- as.integer(nyet)
  y[nyet] <- NA
  return(y)
}

# ifelse(test, -x, x)
ifelseNegPos <- function(test, x) {
  y <- x
  y[test] <- -x[test]
  return(y)
}

# ..................

"%orifnull%" <- function(a, b) {
  if(!is.null(a)) return(a)
  # b is evaluated only now
  return(b)
}

blockdiagmatrix <- function(...) {
  x <- list(...)
  if(!all(unlist(lapply(x, is.matrix))))
    stop("Some of the arguments are not matrices", call.=FALSE)
  nr <- unlist(lapply(x, nrow))
  nc <- unlist(lapply(x, ncol))
  result <- matrix(0, sum(nr), sum(nc))
  rownames(result) <- unlist(lapply(x, rownames))
  colnames(result) <- unlist(lapply(x, colnames))
  rowend <- cumsum(nr)
  rowstart <- c(0, rowend) + 1
  colend <- cumsum(nc)
  colstart <- c(0, colend) + 1
  for(i in seq_along(x))
    result[ (rowstart[i]):(rowend[i]) , (colstart[i]):(colend[i])] <- x[[i]]
  return(result)
}

blockdiagarray <- function(...) {
  x <- list(...)
  if(!all(unlist(lapply(x, is.array))))
    stop("Some of the arguments are not arrays", call.=FALSE)
  dims <- lapply(x, dim)
  dims1 <- unlist(lapply(dims, "[", i=1))
  if(length(dim1 <- unique(dims1)) > 1)
    stop("Arrays have different extents in first dimension")
  dims2 <- unlist(lapply(dims, "[", i=2))
  dims3 <- unlist(lapply(dims, "[", i=3))
  result <- array(0, dim=c(dim1, sum(dims2), sum(dims3)))
  dn <- lapply(x, dimnames)
  dimnames(result)[[2]] <- unlist(lapply(dn, "[[", i=2))
  dimnames(result)[[3]] <- unlist(lapply(dn, "[[", i=3))
  rowend <- cumsum(dims2)
  rowstart <- c(0, rowend) + 1
  colend <- cumsum(dims3)
  colstart <- c(0, colend) + 1
  for(i in seq_along(x))
    result[ , (rowstart[i]):(rowend[i]) , (colstart[i]):(colend[i])] <- x[[i]]
  return(result)
}

lty2char <- function(i) {
  if(is.numeric(i)) c("blank", "solid", "dashed", "dotted",
                      "dotdash", "longdash", "twodash")[(i %% 7) + 1] else i
}

## convert numeric matrix to character, and blank out lower sub-diagonal.
uptrimat <- function(x) {
  stopifnot(is.matrix(x))
  x[] <- as.character(x)
  x[row(x) > col(x)] <- ""
  return(noquote(x))
}

asNumericMatrix <- function(x) {
  ## workaround for strange artefact of as.matrix.data.frame
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  x
}

prepareTitle <- function(main) {
  ## Count the number of lines in a main title
  ## Convert title to a form usable by plot.owin
  if(is.expression(main)) {
    nlines <- 1
  } else {
    main <- paste(main)
    ## break at newline 
    main <- unlist(strsplit(main, "\n"))
    nlines <- if(sum(nchar(main)) == 0) 0 else length(main)
  }
  return(list(main=main,
              nlines=nlines,
              blank=rep('  ', nlines)))
}

simplenumber <- local({

  iswhole <- function(x, tol=0) { abs(x %% 1) <= tol }

  iszero <- function(x, tol=0) { abs(x) <= tol }
  
  simplenumber <- function(x, unit = "", multiply="*",
                           tol=.Machine$double.eps) {
    ## Try to express x as a simple multiple or fraction
    if(length(x) > 1)
      return(sapply(as.list(x), simplenumber,
                    unit=unit, multiply=multiply, tol=tol))
    s <- if(x < 0) "-" else ""
    x <- abs(x)
    if(unit == "") {
      if(iswhole(x, tol)) return(paste0(s, round(x)))
      for(i in 1:12) {
        if(iswhole(i/x, tol)) return(paste0(s, i, "/", round(i/x)))
        if(iswhole(i*x, tol)) return(paste0(s, round(i*x), "/", i))
      }
    } else {
      if(iszero(x, tol)) return("0")
      if(iszero(x-1, tol)) return(paste0(s,unit))
      if(iswhole(x, tol)) return(paste0(s, round(x), multiply, unit))
      if(iswhole(1/x, tol)) return(paste0(s, unit, "/", round(1/x)))
      for(i in 2:12) {
        if(iswhole(i/x, tol))
          return(paste0(s, i, multiply, unit, "/", round(i/x)))
        if(iswhole(i*x, tol))
          return(paste0(s, round(i*x), multiply, unit, "/", i))
      }
    }
    return(NULL)
  }

  simplenumber
})


fontify <- function(x, font="italic") {
  if(!nzchar(font) || font == "plain")
    return(x)
  if(is.character(x))
    return(paste0(font, "(", x, ")"))
  if(is.expression(x)) {
    if((n <- length(x)) > 0) {
      for(i in 1:n) 
        x[[i]] <- fontify(x[[i]], font)
    }
    return(x)
  }
  if(is.language(x) || is.numeric(x)) 
    return(substitute(f(X), list(f=as.name(font), X=x)))
  if(all(sapply(x, is.language)))
    return(lapply(x, fontify))
  return(NULL)
}

dround <- function(x) {
  round(x, getOption('digits'))
}

there.is.no.try <- function(...) {
  #' do, or do not
  y <- try(..., silent=TRUE)
  if(inherits(y, "try-error")) return(NULL)
  return(y)
}

# require a namespace and optionally check whether it is attached
kraever <- function(package, fatal=TRUE) {
  if(!requireNamespace(package, quietly=TRUE)) {
    if(fatal)
      stop(paste("The package", sQuote(package), "is required"),
           call.=FALSE)
    return(FALSE)
  }
  if(spatstat.options(paste("check", package, "loaded", sep=".")) &&
    !isNamespaceLoaded(package)){
    if(fatal)
      stop(paste("The package", sQuote(package),
                 "must be loaded: please type",
                 sQuote(paste0("library", paren(package)))),
           call.=FALSE)
    return(FALSE)
  }
  return(TRUE)
}

## replace recognise keywords by other keywords
mapstrings <- function(x, map=NULL) {
  if(is.null(map)) return(x)
  x <- as.character(x)
  from <- names(map)
  map <- as.character(map)
  if(sum(nzchar(from)) != length(map))
    stop("input names are missing in map", call.=FALSE)
  if(any(duplicated(from)))
    stop("input names are duplicated in map", call.=FALSE)
  i <- match(x, from)
  hit <- !is.na(i)
  x[hit] <- map[i[hit]]
  return(x)
}

romansort <- local({

  # sort character strings in order of Roman alphabet
  
  romansort <- function(x) {
    if(!is.character(x)) return(sort(x))
    x <- as.vector(x)
    ## convert each 'word' to a vector of single characters
    cc <- strsplit(x, "")
    ## find position of each character in Roman alphabet
    mm <- lapply(cc, match, table=c(letters, LETTERS))
    mmax <- max(unlist(mm), na.rm=TRUE)
    ## encode
    nn <- sapply(mm, powercode, base=mmax)
    ## find ordering
    oo <- order(nn, na.last=TRUE)
    return(x[oo])
  }

  powercode <- function(x, base) sum(x * base^rev((seq_len(length(x))-1)))

  romansort
})

checkbigmatrix <- function(n, m, fatal=FALSE, silent=FALSE) {
  if(n * m <= spatstat.options("maxmatrix"))
    return(TRUE)
  whinge <- paste("Attempted to create binary mask with",
                  n, "*", m, "=", n * m, "entries")
  if(fatal) stop(whinge, call.=FALSE)
  if(!silent) warning(whinge, call.=FALSE)
  return(FALSE)
}

insertinlist <- function(x, i, y) {
  ## insert a possibly longer or shorter list 'y'
  ## into serial position 'i' in list 'x'
  n <- length(x)
  if(n == 0) return(y)
  m <- seq_len(n)
  names(m) <- names(x)
  i <- m[[i]] # convert 'i' to integer index
  stopifnot(length(i) == 1)
  if(n == 1) return(y)
  xleft <- x[seq_len(i-1)]
  xright <- x[i + seq_len(n-i)]
  z <- c(xleft, y, xright)
  return(z)
}

  
