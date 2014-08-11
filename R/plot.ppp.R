#
#	plot.ppp.R
#
#	$Revision: 1.67 $	$Date: 2014/04/26 01:01:31 $
#
#
#--------------------------------------------------------------------------

plot.ppp <- local({

  ## determine symbol map for marks of points
  default.symap.points <- function(x, ..., 
                                  chars=NULL, cols=NULL, 
                                  maxsize=NULL, markscale=NULL) {
    marx <- marks(x)
    if(is.null(marx)) {
      ## null or constant map
      return(symbolmap(..., chars=chars, cols=cols))
    }
    if(!is.null(dim(marx)))
      stop("Internal error: multivariate marks in default.symap.points")

    argnames <- names(list(...))
    shapegiven <- "shape" %in% argnames
    chargiven <- (!is.null(chars)) || ("pch" %in% argnames)
    assumecircles <- !(shapegiven || chargiven)
    sizegiven <- ("size" %in% argnames) ||
                 (("cex" %in% argnames) && !shapegiven)
    
    if(inherits(marx, c("Date", "POSIXt"))) {
      ## ......... marks are dates or date/times .....................
      timerange <- range(marx, na.rm=TRUE)
      if(sizegiven) {
        g <- do.call("symbolmap",
          resolve.defaults(list(range=timerange),
                           list(...),
                           if(assumecircles) list(shape="circles") else list(),
                           list(chars=chars, cols=cols)))
        return(g)
      }
      ## attempt to determine a scale for the marks 
      y <- scaletointerval(marx, 0, 1, timerange)
      y <- y[is.finite(y)]
      if(length(y) == 0) return(symbolmap(..., chars=chars, cols=cols))
      scal <- mark.scale.default(y, as.owin(x), 
                                 markscale=markscale, maxsize=maxsize)
      if(is.na(scal)) return(symbolmap(..., chars=chars, cols=cols))
      ## scale determined
      sizefun <- function(x) { (scal/2) * scaletointerval(x, 0, 1, timerange)} 
      g <- do.call("symbolmap",
                   resolve.defaults(list(range=timerange),
                                    list(...),
                                    list(shape="circles",
                                         size=sizefun)))
      return(g)
    }
    if(is.numeric(marx)) {
      ## ............. marks are numeric values ...................
      marx <- marx[is.finite(marx)]
      if(length(marx) == 0)
        return(symbolmap(..., chars=chars, cols=cols))
      markrange <- range(marx)
      ## 
      if(sizegiven) {
        g <- do.call("symbolmap",
          resolve.defaults(list(range=markrange),
                           list(...),
                           if(assumecircles) list(shape="circles") else list(),
                           list(chars=chars, cols=cols)))
        return(g)
      }
      ## attempt to determine a scale for the marks 
      if(all(markrange == 0))
        return(symbolmap(..., chars=chars, cols=cols))
      scal <- mark.scale.default(marx, as.owin(x), 
                                 markscale=markscale, maxsize=maxsize)
      if(is.na(scal)) return(symbolmap(..., chars=chars, cols=cols))
      ## scale determined
      if(markrange[1] >= 0) {
        ## all marks are nonnegative
        g <- symbolmap(range=markrange,
                         ...,
                         shape="circles",
                         size=function(x) { scal * x/2 },
                         cols=cols)
      } else {
        ## some marks are negative
        g <- symbolmap(range=markrange,
                         ...,
                         shape=function(x) ifelse(x >= 0, "circles", "squares"),
                         size=function(x) { scal * ifelse(x >= 0, x/2, -x) },
                         cols=cols)
      }
      return(g)
    }
    ##  ...........  non-numeric marks .........................
    um <- if(is.factor(marx)) levels(marx) else sort(unique(marx))
    ntypes <- length(um)
    ## resolve parameters 'chars' and 'cols'
    chars <- default.charmap(ntypes, chars)
    if(!is.null(cols))
      cols <- rep.int(cols, ntypes)[1:ntypes]
    g <- symbolmap(inputs=um, ..., chars=chars, cols=cols)
    return(g)
  }
                                  
  default.charmap <- function(n, ch=NULL) {
    if(!is.null(ch))
      return(rep.int(ch, n)[1:n])
    if(n <= 25)
      return(1:n)
    ltr <- c(letters, LETTERS)
    if(n <= 52)
      return(ltr[1:n])
    ## wrapped sequence of letters
    warning("Too many types to display every type as a different character")
    return(ltr[1 + (0:(n - 1) %% 52)])
  }

  ## main function
  plot.ppp <-
    function(x, main, ..., clipwin=NULL,
             chars=NULL, cols=NULL, use.marks=TRUE,
             which.marks=NULL, add=FALSE, type=c("p", "n"), 
             legend=TRUE, leg.side=c("left", "bottom", "top", "right"),
             leg.args=list(),
             symap=NULL, maxsize=NULL, markscale=NULL, zap=0.01, 
             show.window=show.all, show.all=!add, do.plot=TRUE,
             multiplot=TRUE)
{
  if(missing(main))
    main <- short.deparse(substitute(x))

  type <- match.arg(type)

  if(!missing(maxsize) || !missing(markscale))
    warn.once("circlescale",
              "Interpretation of arguments maxsize and markscale",
              "has changed (in spatstat version 1.37-0 and later).",
              "Size of a circle is now measured by its diameter.")
  
  if(!is.null(clipwin))
    x <- x[clipwin]
  
  ## sensible default position
  if(legend) {
    leg.side <- match.arg(leg.side)
    vertical <- (leg.side %in% c("left", "right"))
  }
  
  if(type == "n" || npoints(x) == 0) {
    ## plot the window only
    xwindow <- x$window
    if(do.plot) 
      do.call("plot.owin",
              resolve.defaults(list(xwindow),
                               list(...),
                               list(main=main, invert=TRUE, add=add,
                                    type=if(show.window) "w" else "n")))
    if(is.null(symap)) symap <- symbolmap()
    attr(symap, "bbox") <- as.rectangle(xwindow)
    return(invisible(symap))
  }

  ## ................................................................
  ## Handle multiple columns of marks as separate plots
  ##  (unless add=TRUE or which.marks selects a single column
  ##   or multipage = FALSE)
  if(use.marks && is.data.frame(mx <- marks(x))) {
    implied.all <- is.null(which.marks)
    want.several <- implied.all || is.data.frame(mx <- mx[,which.marks])
    do.several <- want.several && !add && multiplot
    if(do.several) {
      ## generate one plot for each column of marks
      y <- as.listof(lapply(mx, function(z, P) setmarks(P,z), P=x))
      out <- do.call("plot",
                     resolve.defaults(list(x=y, main=main,
                                           show.window=show.window,
                                           do.plot=do.plot),
                                      list(...),
                                      list(equal.scales=TRUE), 
                                      list(legend=legend,
                                           leg.side=leg.side,
                                           leg.args=leg.args),
                                      list(chars=chars, cols=cols,
                                           maxsize=maxsize, markscale=markscale,
                                           zap=zap)))
      return(invisible(out))
    } 
    if(is.null(which.marks)) {
      which.marks <- 1
      if(do.plot) message("Plotting the first column of marks")
    }
  }
  
  ## ............... unmarked, or single column of marks ....................

  ## Determine symbol map and mark values to be used
  y <- x
  if(!is.marked(x) || !use.marks) {
    ## Marks are not mapped.
    marx <- NULL
    if(is.null(symap)) symap <- symbolmap(..., chars=chars, cols=cols)
  } else {
    ## Marked point pattern
    marx <- marks(y, dfok=TRUE)
    if(is.data.frame(marx)) {
      ## select column or take first colum
      marx <- marx[, which.marks]
      y <- setmarks(y, marx)
    }
    ok <- complete.cases(as.data.frame(x))
    if(!any(ok)) {
      warning("All mark values are NA; plotting locations only.")
      if(is.null(symap)) symap <- symbolmap()
    } else if(any(!ok)) {
      warning(paste("Some marks are NA;",
                    "corresponding points are omitted."))
      x <- x[ok]
      y <- y[ok]
      marx <- marks(y)
    }
    ## apply default symbol map
    if(is.null(symap))
      symap <- default.symap.points(y, chars=chars, cols=cols, 
                                  maxsize=maxsize, markscale=markscale,
                                  ...)
  }
  gtype <- symbolmaptype(symap)

  ## Determine bounding box for main plot
  BB <- as.rectangle(x)
  sick <- inherits(x, "ppp") && !is.null(rejects <- attr(x, "rejects"))
  if(sick) {
    ## Get relevant parameters
    par.direct <- list(main=main, use.marks=use.marks,
                   maxsize=maxsize, markscale=markscale)
    par.rejects <- resolve.1.default(list(par.rejects=list(pch="+")),
                                     list(...))
    par.all <- resolve.defaults(par.rejects, par.direct)
    rw <- resolve.defaults(list(...), list(rejectwindow=NULL))$rejectwindow
    ## determine window for rejects
    rwin <-
      if(is.null(rw))
        rejects$window
      else if(is.logical(rw) && rw)
        rejects$window
      else if(inherits(rw, "owin"))
        rw
      else if(is.character(rw)) {
        switch(rw,
               box={boundingbox(rejects, x)},
               ripras={ripras(c(rejects$x, x$x), c(rejects$y, x$y))},
               stop(paste("Unrecognised option: rejectwindow=", rw)))
      } else stop("Unrecognised format for rejectwindow")
    if(is.null(rwin))
      stop("Selected window for rejects pattern is NULL")
    BB <- boundingbox(BB, as.rectangle(rwin))
  }

  ## Augment bounding box with space for legend, if appropriate
  legend <- legend && (symbolmaptype(symap) != "constant") 
  if(legend) {
    ## guess maximum size of symbols
    maxsize <- invoke.symbolmap(symap, marx,
                                  corners(as.rectangle(x)),
                                  add=add, do.plot=FALSE)
    leg.args <- append(list(side=leg.side, vertical=vertical), leg.args)
    ## draw up layout
    legbox <- do.call.matched(plan.legend.layout,
                              append(list(B=BB, size = 1.5 * maxsize,
                                          started=FALSE, map=symap),
                                     leg.args))
    ## bounding box for everything
    BB <- legbox$A
  }

  ## return now if not plotting
  attr(symap, "bbox") <- BB
  if(!do.plot)
    return(invisible(symap))
    
  ## ............. start plotting .......................
  plot(BB, type="n", add=add, main="           ", show.all=show.all)

  if(sick) {
    if(show.window) {
      ## plot windows
      if(!is.null(rw)) {
        ## plot window for rejects
        rwinpardefault <- list(lty=2,lwd=1,border=1)
        rwinpars <-
          resolve.defaults(par.rejects, rwinpardefault)[names(rwinpardefault)]
        do.call("plot.owin", append(list(rwin, add=TRUE), rwinpars))
      }
      ## plot window of main pattern
      do.call("plot.owin",
              resolve.defaults(list(x$window, add=TRUE),
                               list(...),
                               list(invert=TRUE)))
    }
    ## plot reject points
    do.call("plot.ppp", append(list(rejects, add=TRUE), par.all))
    warning(paste(rejects$n, "illegal points also plotted"))
    ## the rest is added
    add <- TRUE
  }

  ## Now convert to bona fide point pattern
  x <- as.ppp(x)
  xwindow <- x$window

  ## Plot observation window
  if(show.window) {
    do.call("plot.owin",
            resolve.defaults(list(x=xwindow, add=TRUE),
                           list(...),
                           list(invert=TRUE,
                                main=main, show.all=show.all)))
  } else if(show.all) fakemaintitle(as.rectangle(xwindow), main, ...)
  
  ## plot symbols ##
  invoke.symbolmap(symap, marx, x, add=TRUE)

  ## add legend
  if(legend) {
    b <- legbox$b
    do.call("plot",
            append(list(x=symap, main="", add=TRUE,
                        xlim=b$xrange, ylim=b$yrange),
                   leg.args))
  }
  
  return(invisible(symap))
}

plot.ppp

})


mark.scale.default <- function(marx, w, markscale=NULL, maxsize=NULL) {
  ## establish values of markscale, maxsize
  if(!is.null(maxsize) && !is.null(markscale))
    stop("Only one of maxsize and markscale should be given")
  if(is.null(maxsize) && is.null(markscale)) {
    ## if BOTH are absent, enforce the spatstat defaults
    ## (which could also be null)
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
    ## guess appropriate max physical size of symbols
    bb <- as.rectangle(w)
    maxsize <- 1.4/sqrt(pi * length(marx)/area.owin(bb))
    maxsize <- min(maxsize, diameter(bb) * 0.07)
    ## updated: maxsize now represents *diameter*
    maxsize <- 2 * maxsize
  } else stopifnot(maxsize > 0)
  
  # Examine mark values
  maxabs <- max(abs(marx))
  tiny <- (maxabs < 4 * .Machine$double.eps)
  if(tiny)
    return(NA)
  else 
   return(maxsize/maxabs)
}

fakemaintitle <- function(bb, main, ...) {
  ## Try to imitate effect of 'title(main=main)' above a specified box
  if(!any(nzchar(main))) return(invisible(NULL))
  bb <- as.rectangle(bb)
  x0 <- mean(bb$xrange)
  y0 <- bb$yrange[2] + length(main) * diff(bb$yrange)/12
  parnames <- c('cex.main', 'col.main', 'font.main')
  parlist <- par(parnames)
  parlist <- resolve.defaults(list(...), parlist)[parnames]
  names(parlist) <- c('cex', 'col', 'font')
  do.call.matched("text.default",
                  resolve.defaults(list(x=x0, y=y0, labels=main),
                                   parlist,    list(...)))
  return(invisible(NULL))
}
