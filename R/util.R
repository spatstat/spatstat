#
#    util.S    miscellaneous utilities
#
#    $Revision: 1.112 $    $Date: 2013/04/25 06:37:43 $
#
#  (a) for matrices only:
#
#    matrowany(X) is equivalent to apply(X, 1, any)
#    matrowall(X) "   "  " "  "  " apply(X, 1, all)
#    matcolany(X) "   "  " "  "  " apply(X, 2, any)
#    matcolall(X) "   "  " "  "  " apply(X, 2, all)
#
#  (b) for 3D arrays only:
#    apply23sum(X)  "  "   "  " apply(X, c(2,3), sum)
#
#  (c) weighted histogram
#    whist()
#
#  (d) for matrices:
#    matrixsample()
#           subsamples or supersamples a matrix
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
      cell <- factor(cell, levels=0:N)
      # compute weighted histogram
      if(is.null(weights))
        h <- table(cell)
      else 
        h <- unlist(lapply(split(weights, cell), sum, na.rm=TRUE))
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
  xx <- raster.x(masque)
  yy <- raster.y(masque)
  xx <- xx[masque$m]
  yy <- yy[masque$m]
  return(ppp(xx, yy, W))
}

# text magic

commasep <- function(x, join="and") {
  px <- paste(x)
  nx <- length(px)
  if(nx <= 1) return(px)
  commas <- c(rep(", ", length(px)-2),
              paste("", join, ""),
              "")
  return(paste(paste(px, commas, sep=""), collapse=""))
}

paren <- function(x, type="(") {
  switch(type,
         "(" = {
           out <- paste("(", x, ")", sep="")
         },
         "[" = {
           out <- paste("[", x, "]", sep="")
         },
         "{" = {
           out <- paste("{", x, "}", sep="")
         })
  out
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

# equivalent to rev(cumsum(rev(x)))

revcumsum <- function(x) {
  n <- length(x)
  if(identical(storage.mode(x), "integer")) {
    z <- .C("irevcumsum",
            x=as.integer(x),
            as.integer(n),
            PACKAGE="spatstat")
    return(z$x)
  } else {
    z <- .C("drevcumsum",
            x=as.double(x),
            as.integer(n),
            PACKAGE="spatstat")
    return(z$x)
  }
}

prolongseq <- function(x, newrange) {
  stopifnot(length(newrange) == 2 && newrange[1] < newrange[2])
  stopifnot(length(x) >= 2)
  dx <- diff(x)
  if(any(dx <= 0))
    stop("x must be an increasing sequence")
  if(diff(range(dx)) > 0.01 * abs(mean(dx)))
    stop("x must be evenly spaced")
  dx <- mean(dx)

  # add or trim data to left
  if(x[1] > newrange[1]) {
    leftbit <- seq(from=x[1], to=newrange[1], by= -dx) 
    x <- c(rev(leftbit), x[-1])
  } else 
    x <- x[x >= newrange[1]]

  # add or trim data to right
  nx <- length(x)
  if(newrange[2] > x[nx]) {
    rightbit <- seq(from=x[nx], to=newrange[2], by= dx)
    x <- c(x[-nx], rightbit)
  } else 
    x <- x[x <= newrange[2]]

  return(x)
}

intersect.ranges <- function(a, b, fatal=TRUE) {
  lo <- max(a[1],b[1])
  hi <- min(a[2],b[2])
  if(lo >= hi) {
    if(fatal) stop("Intersection is empty")
    else return(NULL)
  }
  return(c(lo, hi))
}

inside.range <- function(x, r) {
  stopifnot(length(r) == 2 && r[1] < r[2])
  return(x >= r[1] & x <= r[2])
}

prettyinside <- function(x, ...) {
  r <- range(x, na.rm=TRUE)
  p <- pretty(x, ...)
  ok <- inside.range(p, r)
  return(p[ok])
}

check.range <- function(x, fatal=TRUE) {
  xname <- deparse(substitute(x))
  if(identical(x, range(x, na.rm=TRUE)))
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

assign(".Spatstat.ProgressBar", NULL, envir = .spEnv)
assign(".Spatstat.ProgressData", NULL, envir = .spEnv)

progressreport <- function(i, n, every=min(100,max(1, ceiling(n/100))),
                           nperline=min(charsperline,
                             every * ceiling(charsperline /(every+3))),
                           charsperline=60,
                           style=spatstat.options("progress")) {
  missevery <- missing(every)
  if(i > n) {
    warning(paste("progressreport called with i =", i, "> n =", n))
    return(invisible(NULL))
  }
  switch(style,
         txtbar={
           if(i == 1) {
             # initialise text bar
             assign(".Spatstat.ProgressBar",
                    txtProgressBar(1, n, 1, style=3),
                    envir = .spEnv)
           } else {
             # get text bar
             pbar <- get(".Spatstat.ProgressBar", envir = .spEnv)
             # update 
             setTxtProgressBar(pbar, i)
             if(i == n) {
               close(pbar)
               assign(".Spatstat.ProgressBar", NULL, envir = .spEnv)
             } 
           }
         },
         tty={
           now <- proc.time()
           if(i == 1) {
             # Initialise stuff
             if(missevery && every > 1 && n > 10) {
               every <- niceround(every)
               nperline <- min(charsperline,
                               every * ceiling(charsperline /(every+3)))
             }
             showtime <- FALSE
             showevery <- n
             assign(".Spatstat.ProgressData",
                    list(every=every, nperline=nperline,
                         starttime=now,
                         showtime=FALSE, showevery=n),
                    envir=.spEnv)
           } else {
             pd <- get(".Spatstat.ProgressData", envir=.spEnv)
             if(is.null(pd))
               stop(paste("progressreport called with i =", i, "before i = 1"))
             every     <- pd$every
             nperline  <- pd$nperline
             showtime  <- pd$showtime
             showevery <- pd$showevery
             if(i < n) {
               # estimate time remaining
               starttime <- pd$starttime
               elapsed <- now - starttime
               elapsed <- unname(elapsed[3])
               rate <- elapsed/(i-1)
               remaining <- rate * (n-i)
               if(!showtime) {
                 # show time remaining if..
                 if(rate > 20) {
                   # .. rate is very slow
                   showtime <- TRUE
                   showevery <- 1
                 } else if(remaining > 180) {
                   # ... more than 3 minutes remaining
                   showtime <- TRUE
                   showevery <- every
                   aminute <- ceiling(60/rate)
                   if(aminute < showevery) 
                     showevery <- min(niceround(aminute), showevery)
                 }
               }
               assign(".Spatstat.ProgressData",
                    list(every=every, nperline=nperline,
                         starttime=starttime,
                         showtime=showtime, showevery=showevery),
                    envir=.spEnv)
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
  return(invisible(NULL))
}
  
numalign <- function(i, nmax, zero="0") {
  stopifnot(i <= nmax)
  nplaces <- as.integer(ceiling(log10(nmax+1)))
  out <- blank <- paste(rep(zero, nplaces), collapse="")
  istring <- paste(i)
  ilen <- nchar(istring)
  substr(out, nplaces-ilen+1, nplaces) <- istring
  return(out)
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

check.nvector <- function(v, npoints, fatal=TRUE, things="data points",
                          naok=FALSE) {
  # vector of numeric values for each point/thing
  vname <- sQuote(deparse(substitute(v)))
  whinge <- NULL
  if(!is.numeric(v))
    whinge <- paste(vname, "is not numeric")
  else if(!is.atomic(v) || !is.null(dim(v)))  # vector with attributes
    whinge <- paste(vname, "is not a vector")
  else if(length(v) != npoints)
    whinge <- paste("The length of", vname,
                    "should equal the number of", things)
  else if(!naok && any(is.na(v)))
    whinge <- paste("Some values of", vname, "are NA or NaN")
  #
  if(!is.null(whinge)) {
    if(fatal) stop(whinge) else return(FALSE)
  }
  return(TRUE)
}

check.nmatrix <- function(m, npoints, fatal=TRUE, things="data points",
                          naok=FALSE, squarematrix=TRUE) {
  # matrix of values for each thing or each pair of things
  mname <- sQuote(deparse(substitute(m)))
  whinge <- NULL
  if(!is.matrix(m))
    whinge <- paste(mname, "should be a matrix")
  else if(squarematrix && (nrow(m) != ncol(m)))
    whinge <- paste(mname, "should be a square matrix")
  else if(!naok && any(is.na(m)))
    whinge <- paste("Some values of", mname, "are NA or NaN")
  else if(nrow(m) != npoints)
    whinge <- paste("Number of rows in", mname,
               "does not match number of", things)
  #
  if(!is.null(whinge)) {
    if(fatal) stop(whinge) else return(FALSE)
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

firstfactor <- function(x) {
  stopifnot(is.data.frame(x) || is.hyperframe(x))
  isfac <- unlist(lapply(as.list(x), is.factor))
  if(!any(isfac)) 
    return(NULL)
  return(x[, min(which(isfac)), drop=TRUE])
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
  ifelse(is.finite(x) & x != 0, x * a, x)
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

short.deparse <- function(x, maxlen=60) {
  deparse(x,
          nlines=1,
          width.cutoff=maxlen,
          control="delayPromises")
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
  if(any(duplicated(result)))
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
  codetime <- function(x, hms=TRUE) {
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

# print names and version numbers of libraries loaded

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

dropifsingle <- function(x) if(length(list) == 1) x[[1]] else x

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
  if(is.numeric(x)) print(as.numeric(x), ...) else NextMethod("print")
  taken <- summary(attr(x, "timetaken"))[[1]]
  cat(paste("\nTime taken:", codetime(taken), "\n"))
  return(invisible(NULL))
}
