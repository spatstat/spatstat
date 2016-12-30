#
#    util.R    miscellaneous utilities
#
#    $Revision: 1.236 $    $Date: 2016/12/30 03:24:37 $
#

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

onecolumn <- function(m) {
  switch(markformat(m),
         none=stop("No marks provided"),
         vector=m,
         dataframe=m[,1, drop=TRUE],
         NA)
}


checkbigmatrix <- function(n, m, fatal=FALSE, silent=FALSE) {
  if(n * m <= spatstat.options("maxmatrix"))
    return(TRUE)
  whinge <- paste("Attempted to create binary mask with",
                  n, "*", m, "=", n * m, "entries")
  if(fatal) stop(whinge, call.=FALSE)
  if(!silent) warning(whinge, call.=FALSE)
  return(FALSE)
}


## ........... progress reports .....................

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

  IterationsPerLine <- function(charsperline, n, every, tick,
                                showtime, showevery) {
    # Calculate number of iterations that triggers a newline.
    # A dot is printed every 'tick' iterations
    # Iteration number is printed every 'every' iterations.
    # If showtime=TRUE, the time is shown every 'showevery' iterations
    # where showevery \in {1, every, n}.
    chars.report <- max(1, ceiling(log10(n)))
    if(showtime) {
      chars.time <- nchar(' [etd 12:00:00] ')
      timesperreport <- if(showevery == 1) every else
                        if(showevery == every) 1 else 0
      chars.report <- chars.report + timesperreport * chars.time
    }
    chars.ticks <- floor((every-1)/tick)
    chars.block <- chars.report + chars.ticks
    nblocks <- max(1, floor(charsperline/chars.block))
    nperline <- nblocks * every
    leftover <- charsperline - nblocks * chars.block
    if(leftover > 0)
      nperline <- nperline + min(leftover * tick, showevery - 1)
    return(nperline)
  }
  
  progressreport <- function(i, n,
                             every=min(100,max(1, ceiling(n/100))),
                             tick=1,
                             nperline=NULL,
                             charsperline=getOption("width"),
                             style=spatstat.options("progress"),
                             showtime=NULL,
                             state=NULL) {
    missevery <- missing(every)
    nperline.fixed <- !is.null(nperline)
    showtime.optional <- is.null(showtime)
    if(showtime.optional) showtime <- FALSE # initialise only
    if(i > n) {
      warning(paste("progressreport called with i =", i, "> n =", n))
      return(invisible(NULL))
    }
    if(style == "tk" && !requireNamespace("tcltk")) {
      warning("tcltk is unavailable; switching to style='txtbar'", call.=FALSE)
      style <- "txtbar"
    }
    if(is.null(state) && style != "tty")
      stop(paste("Argument 'state' is required when style =",sQuote(style)),
           call.=FALSE)
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
             if(i == 1 || is.null(state)) {
               ## Initialise stuff
               if(missevery && every > 1 && n > 10) 
                 every <- niceround(every)
               showevery <- if(showtime) every else n
               if(!nperline.fixed) 
                 nperline <- IterationsPerLine(charsperline, n, every, tick,
                                               showtime, showevery)
               state <- Put("ProgressData",
                            list(every=every,
                                 tick=tick,
                                 nperline=nperline,
                                 starttime=now,
                                 showtime=showtime,
                                 showevery=showevery,
                                 nperline.fixed=nperline.fixed,
                                 showtime.optional=showtime.optional),
                            state)
             } else {
               pd <- Get("ProgressData", state)
               if(is.null(pd))
                 stop(paste("progressreport called with i =", i,
                            "before i = 1"))
               every     <- pd$every
               tick      <- pd$tick
               nperline  <- pd$nperline
               showtime  <- pd$showtime
               showevery <- pd$showevery
               showtime.optional <- pd$showtime.optional
               nperline.fixed    <- pd$nperline.fixed
               if(i < n) {
                 if(showtime || showtime.optional) {
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
                     # update number of iterations per line
                     if(showtime && !nperline.fixed) 
                       nperline <- IterationsPerLine(charsperline,
                                                     n, every, tick,
                                                     showtime, showevery)
                   }
                 }
                 state <- Put("ProgressData",
                              list(every=every,
                                   tick=tick,
                                   nperline=nperline,
                                   starttime=starttime,
                                   showtime=showtime,
                                   showevery=showevery,
                                   nperline.fixed=nperline.fixed,
                                   showtime.optional=showtime.optional),
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
               else if(i %% tick == 0)
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
  


multiply.only.finite.entries <- function(x, a) {
  # In ppm a potential value that is -Inf must remain -Inf
  # and a potential value that is 0 multiplied by NA remains 0
  y <- x
  ok <- is.finite(x) & (x != 0)
  y[ok] <- a * x[ok]
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



# ..................

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

requireversion <- function(pkg, ver) {
  pkgname <- deparse(substitute(pkg))
  v <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), 
                fields="Version")
  if(package_version(v) < ver)
    stop(paste("Package",
               sQuote(pkgname),
               "is out of date: version >=",
               ver,
               "is needed"),
         call.=FALSE)
  invisible(NULL)
}

spatstatDiagnostic <- function(msg) {
  cat("-----------------------------\n")
  cat(paste(" >>> Spatstat Diagnostic: ", msg, "<<<\n"))
  cat("-----------------------------\n")
  invisible(NULL)
}

