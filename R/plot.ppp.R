#
#	plot.ppp.R
#
#	$Revision: 1.52 $	$Date: 2013/12/04 11:10:17 $
#
#
#--------------------------------------------------------------------------

plot.ppp <-
  function(x, main, ..., chars=NULL, cols=NULL, use.marks=TRUE,
           which.marks=NULL, add=FALSE, type=c("p", "n"), 
           maxsize=NULL, markscale=NULL, zap=0.01)
{
  if(missing(main))
    main <- short.deparse(substitute(x))

  type <- match.arg(type)
  
  if(type == "n") {
    # plot the window only
    do.call("plot.owin",
            resolve.defaults(list(x$window),
                             list(...),
                             list(main=main, invert=TRUE, add=add)))
    return(invisible(NULL))
  }
  
# Handle multiple columns of marks as separate plots
#  (unless add=TRUE or which.marks selects a single column)
  if(use.marks && is.data.frame(mx <- marks(x))) {
    implied.all <- is.null(which.marks)
    do.several <- implied.all || is.data.frame(mx <- mx[,which.marks])
    if(add && implied.all) {
      message("Plotting the first column of marks")
      which.marks <- 1
    } else if(!add && do.several) {
      y <- as.listof(lapply(mx, function(z, P) setmarks(P,z), P=x))
      out <- do.call("plot",
                     resolve.defaults(list(x=y, main=main),
                                      list(...),
                                      list(chars=chars, cols=cols,
                                           maxsize=maxsize, markscale=markscale,
                                           zap=zap)))
      if(is.null(out)) return(invisible(NULL)) else return(out)
    } 
  }

# First handle `rejected' points
  sick <- inherits(x, "ppp") && !is.null(rejects <- attr(x, "rejects"))
  if(sick) {
    # get any parameters
    par.direct <- list(main=main, use.marks=use.marks,
                   maxsize=maxsize, markscale=markscale)
    par.rejects.default <- list(pch="+")
    par.rejects <- resolve.defaults(list(...),
                                    list(par.rejects=par.rejects.default))$par.rejects
    par.rejects <- resolve.defaults(par.rejects, par.rejects.default)
    par.all <- resolve.defaults(par.rejects, par.direct)
    rw <- resolve.defaults(list(...), list(rejectwindow=NULL))$rejectwindow
    # determine window for rejects
    rwin <-
      if(is.null(rw))
        rejects$window
      else if(is.logical(rw) && rw)
        rejects$window
      else if(inherits(rw, "owin"))
        rw
      else if(is.character(rw)) {
        switch(rw,
               box={bounding.box(rejects, x)},
               ripras={ripras(c(rejects$x, x$x), c(rejects$y, x$y))},
               stop(paste("Unrecognised option: rejectwindow=", rw)))
      } else stop("Unrecognised format for rejectwindow")
    if(is.null(rwin))
      stop("Selected window for rejects pattern is NULL")
    # Create suitable space
    plot(rejects$window, add=add, type="n", main="")
    if(!add)
      title(main=main)
    # plot rejects window if commanded
    if(!is.null(rw)) {
      rwinpardefault <- list(lty=2,lwd=1,border=1)
      rwinpars <-
        resolve.defaults(par.rejects, rwinpardefault)[names(rwinpardefault)]
      do.call("plot.owin", append(list(rwin, add=TRUE), rwinpars))
    }
    # plot window of main pattern
    do.call("plot.owin",
            resolve.defaults(list(x$window, add=TRUE),
                             list(...),
                             list(invert=TRUE)))
    # plot points
    do.call("plot.ppp", append(list(rejects, add=TRUE), par.all))
    warning(paste(rejects$n, "illegal points also plotted"))
    # the rest is added
    add <- TRUE
  }

# Now convert to bona fide point pattern
  x <- as.ppp(x)
  xwindow <- x$window
  marked <- is.marked(x, dfok=TRUE, na.action="ignore")

# Plot observation window
  if(!add)
    do.call("plot.owin",
            resolve.defaults(list(xwindow),
                             list(...),
                             list(invert=TRUE, main=main)))
    
  if(x$n == 0)
    return(invisible())

# Handle plot parameters
  explicit <- list()
  if(!is.null(cols))
    explicit <- append(explicit, list(cols=cols))
  if(!is.null(chars))
    explicit <- append(explicit, list(chars=chars))
    
  defaults <- spatstat.options("par.points")

# Prepare to plot points
  
  smartpoints <- function(xx, yy, ...,
                          index=1, col=NULL, pch=NULL, cols=NULL, chars=NULL) {
    if(!is.null(cols))
      col <- cols[index]
    if(is.null(pch) && !is.null(chars))
      pch <- chars[index]
    do.call.matched("points",
            resolve.defaults(list(x=list(x=xx, y=yy), ...),
                             if(!is.null(col)) list(col=col) else NULL,
                             if(!is.null(pch)) list(pch=pch) else NULL),
                    extrargs=c("col", "pch", "type", "bg", "cex", "lwd", "lty"))
  }

  if(!marked || !use.marks) {
    do.call("smartpoints",
            resolve.defaults(list(xx=x$x, yy=x$y),
                             explicit,
                             list(...),
                             spatstat.options("par.points")))
    return(invisible())
  }

  # marked point pattern

  marx <- marks(x, dfok=TRUE)

  if(is.data.frame(marx)) {
    # select column or take first colum
    marx <- marx[, which.marks]
  }

  # check there are some valid marks!
  ok <- !is.na(marx)
  if(all(!ok)) {
    warning("All mark values are NA; plotting locations only.")
    do.call("smartpoints",
            resolve.defaults(list(xx=x$x, yy=x$y),
                             explicit,
                             list(...),
                             spatstat.options("par.points")))
    return(invisible())
  }

  # otherwise ignore invalid marks
  if(!all(ok)) {
    warning(paste("Some marks are NA;",
                    "corresponding points are omitted."))
    x <- x[ok]
    marx <- marx[ok]
  }

  
  ################  convert POSIX times to real numbers ###########

  if(marks.are.times <- inherits(marx, "POSIXt")) {
    marx <- as.POSIXct(marx)
    tzone <- attr(marx, "tzone")
    earliest.time <- min(marx)
    marx <- as.numeric(marx - earliest.time)
  }
  
  ################  real-valued marks ############################

  if(is.numeric(marx)) {

    ok <- is.finite(marx)

    if(!all(ok)) {
      warning(paste("Some marks are infinite",
                    "corresponding points are omitted."))
      x <- x[ok]
      marx <- marx[ok]
    }

    scal <- mark.scale.default(marx, xwindow,
                               markscale=markscale, maxsize=maxsize)
    if(is.na(scal)) {
      # data cannot be scaled successfully;
      # plot as points
      do.call("smartpoints",
              resolve.defaults(list(x$x, x$y),
                               explicit,
                               list(...),
                               spatstat.options("par.points")))
      return(invisible())
    }
    # scale determined.
    # Apply the scaling
    ms <- marx * scal 

    # Finally, plot them..
    absmarx <- abs(marx)
    tiny <- (absmarx <= zap * max(absmarx))
    neg <- (marx < 0) & !tiny
    pos <- (marx > 0) & !tiny
    # plot positive values as circles
    if(any(pos)) 
      do.call("symbols",
              resolve.defaults(
                               list(x$x[pos], x$y[pos]),
                               list(circles = ms[pos]),
                               list(inches = FALSE, add = TRUE),
                               if(!is.null(cols)) list(fg=cols[1]) else NULL,
                               list(...)))
    # plot negative values as squares
    if(any(neg))
      do.call("symbols",
              resolve.defaults(
                               list(x$x[neg], x$y[neg]),
                               list(squares = - ms[neg]),
                               list(inches = FALSE, add = TRUE),
                               if(!is.null(cols)) list(fg=cols[1]) else NULL,
                               list(...)))
    # return a plottable scale bar
    mr <- range(marx)
    mp.value <- if(is.na(scal)) mr[1] else pretty(mr)
    mp.plotted <- mp.value * scal
    if(marks.are.times)
      mp.value <- as.POSIXct(mp.value,
                             tz=tzone, origin=earliest.time)
    names(mp.plotted) <- paste(mp.value)
    return(mp.plotted)
  }

  ##################### non-numeric marks ###############################

  um <- if(is.factor(marx))
    levels(marx)
  else
    sort(unique(marx))

  ntypes <- length(um)
  
  if(is.null(chars)) {
    if(ntypes <= 25) {
      # numerical 'pch' 
      chars <- 1:ntypes
    } else {
      # letters
      ltr <- c(letters, LETTERS)
      if(ntypes <= 52) {
        chars <- ltr[1:ntypes]
      } else {
        # wrapped sequence of letters
        warning("There are too many types to display every type as a different character")
        chars <- ltr[1 + (0:(ntypes - 1) %% 52)]
      }
    }
  }
  else if((nchars <- length(chars)) != ntypes) {
    if(nchars != 1)
      stop(paste("length of",
                 sQuote("chars"),
                 "is not equal to the number of types"))
    else
      explicit$chars <- chars <- rep.int(chars, ntypes)
  }

  if(!is.null(cols) && ((ncols <- length(cols)) != ntypes)) {
    if(ncols != 1)
      stop(paste("length of",
                 sQuote("cols"),
                 "is not equal to the number of types"))
    else
      explicit$cols <- cols <- rep.int(cols, ntypes)
  }
    
  for(i in seq_along(um)) {
    relevant <- (marx == um[i])
    if(any(relevant))
      do.call("smartpoints",
              resolve.defaults(list(x$x[relevant], x$y[relevant]),
                               list(pch = chars[i]),
                               explicit,
                               list(index=i),
                               list(...),
                               spatstat.options("par.points")))
  }
  names(chars) <- um
  if(length(chars) < 20)
    return(chars)
  else
    return(invisible(chars))
}


mark.scale.default <- function(marx, w, markscale=NULL, maxsize=NULL) {
    # establish values of markscale, maxsize
  if(!is.null(maxsize) && !is.null(markscale))
    stop("Only one of maxsize and markscale should be given")
  if(is.null(maxsize) && is.null(markscale)) {
    # if BOTH are absent, enforce the spatstat defaults
    # (which could also be null)
    pop <- spatstat.options("par.points")
    markscale <- pop$markscale
    maxsize   <- pop$maxsize
  }
  # Now check whether markscale is fixed
  if(!is.null(markscale)) {
    stopifnot(markscale > 0)
    return(markscale)
  }
  # Usual case: markscale is to be determined from maximum physical size
  if(is.null(maxsize)) {
      # guess appropriate max physical size of symbols
    bb <- as.rectangle(w)
    maxsize <- 1.4/sqrt(pi * length(marx)/area.owin(bb))
    maxsize <- min(maxsize, diameter(bb) * 0.07)
  } else stopifnot(maxsize > 0)
  
  # Examine mark values
  maxabs <- max(abs(marx))
  tiny <- (maxabs < 4 * .Machine$double.eps)
  if(tiny)
    return(NA)
  else 
    return(maxsize/maxabs)
}
