#
#	plot.ppp.R
#
#	$Revision: 1.87 $	$Date: 2016/06/14 02:26:13 $
#
#
#--------------------------------------------------------------------------

plot.ppp <- local({

  transparencyfun <- function(n) {
    if(n <= 100) 1 else (0.2 + 0.8 * exp(-(n-100)/1000))
  }
  
  ## determine symbol map for marks of points
  default.symap.points <- function(x, ..., 
                                  chars=NULL, cols=NULL, 
                                  maxsize=NULL, meansize=NULL, markscale=NULL) {
    marx <- marks(x)
    if(is.null(marx)) {
      ## null or constant symbol map
      ## consider using transparent colours
      if(is.null(cols) &&
         !any(c("col", "fg", "bg") %in% names(list(...))) &&
         (nx <- npoints(x)) > 100 &&
         identical(dev.capabilities()$semiTransparency, TRUE) &&
         spatstat.options("transparent"))
        cols <- rgb(0,0,0,transparencyfun(nx))
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
      shapedefault <- if(!assumecircles) list() else list(shape="circles")
      if(sizegiven) {
        g <- do.call(symbolmap,
          resolve.defaults(list(range=timerange),
                           list(...),
                           shapedefault,
                           list(chars=chars, cols=cols)))
        return(g)
      }
      ## attempt to determine a scale for the marks 
      y <- scaletointerval(marx, 0, 1, timerange)
      y <- y[is.finite(y)]
      if(length(y) == 0) return(symbolmap(..., chars=chars, cols=cols))
      scal <- mark.scale.default(y, as.owin(x), 
                                 markscale=markscale, maxsize=maxsize,
                                 meansize=meansize, 
                                 characters=chargiven)
      if(is.na(scal)) return(symbolmap(..., chars=chars, cols=cols))
      ## scale determined
      sizefun <- function(x, scal=1) {
        (scal/2) * scaletointerval(x, 0, 1, timerange)
      }
      formals(sizefun)[[2]] <- scal  ## ensures value of 'scal' is printed
      ##
      g <- do.call(symbolmap,
                   resolve.defaults(list(range=timerange),
                                    list(...),
                                    shapedefault,
                                    list(size=sizefun)))
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
        g <- do.call(symbolmap,
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
                                 markscale=markscale, maxsize=maxsize,
                                 meansize=meansize,
                                 characters=chargiven)
      if(is.na(scal)) return(symbolmap(..., chars=chars, cols=cols))
      ## scale determined
      if(markrange[1] >= 0) {
        ## all marks are nonnegative
        shapedefault <-
          if(!assumecircles) list() else list(shape="circles")
        cexfun <- function(x, scal=1) { scal * x }
        circfun <- function(x, scal=1) { scal * x/2 }
        formals(cexfun)[[2]] <- formals(circfun)[[2]] <- scal
        sizedefault <-
          if(sizegiven) list() else
          if(chargiven) list(cex=cexfun) else list(size=circfun)
      } else {
        ## some marks are negative
        shapedefault <-
          if(!assumecircles) list() else
          list(shape=function(x) { ifelse(x >= 0, "circles", "squares") })
        cexfun <- function(x, scal=1) { scal * abs(x) }
        circfun <- function(x, scal=1) { scal * ifelse(x >= 0, x/2, -x) }
        formals(cexfun)[[2]] <- formals(circfun)[[2]] <- scal
        sizedefault <-
          if(sizegiven) list() else
          if(chargiven) list(cex=cexfun) else list(size=circfun)
      }
      g <- do.call(symbolmap,
                   resolve.defaults(list(range=markrange),
                                    list(...),
                                    shapedefault,
                                    sizedefault,
                                    list(chars=chars, cols=cols)))
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
             symap=NULL, maxsize=NULL, meansize=NULL, markscale=NULL, zap=0.01, 
             show.window=show.all, show.all=!add, do.plot=TRUE,
             multiplot=TRUE)
{
  if(missing(main))
    main <- short.deparse(substitute(x))

  type <- match.arg(type)

  if(!missing(maxsize) || !missing(markscale) || !missing(meansize))
    warn.once("circlescale",
              "Interpretation of arguments maxsize and markscale",
              "has changed (in spatstat version 1.37-0 and later).",
              "Size of a circle is now measured by its diameter.")

  if(clipped <- !is.null(clipwin)) {
    stopifnot(is.owin(clipwin))
    W <- Window(x)
    clippy <- if(is.mask(W)) intersect.owin(W, clipwin) else edges(W)[clipwin]
    x <- x[clipwin]
  } else clippy <- NULL
  
  ## sensible default position
  legend <- legend && show.all
  if(legend) {
    leg.side <- match.arg(leg.side)
    vertical <- (leg.side %in% c("left", "right"))
  }
  
#  if(type == "n" || npoints(x) == 0) {
#    ## plot the window only
#    xwindow <- x$window
#    if(do.plot) 
#      do.call(plot.owin,
#              resolve.defaults(list(xwindow),
#                               list(...),
#                               list(main=main, invert=TRUE, add=add,
#                                    type=if(show.window) "w" else "n")))
#    if(is.null(symap)) symap <- symbolmap()
#    attr(symap, "bbox") <- as.rectangle(xwindow)
#    return(invisible(symap))
#  }

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
      y <- solapply(mx, setmarks, x=x)
      out <- do.call(plot,
                     resolve.defaults(list(x=y, main=main,
                                           show.window=show.window && !clipped,
                                           do.plot=do.plot,
                                           type=type),
                                      list(...),
                                      list(equal.scales=TRUE),
                                      list(panel.end=clippy),
                                      list(legend=legend,
                                           leg.side=leg.side,
                                           leg.args=leg.args),
                                      list(chars=chars, cols=cols,
                                           maxsize=maxsize,
                                           meansize=meansize,
                                           markscale=markscale,
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
  if(!is.marked(x, na.action="ignore") || !use.marks) {
    ## Marks are not mapped.
    marx <- NULL
    if(is.null(symap))
      symap <- default.symap.points(unmark(x), ..., chars=chars, cols=cols)
  } else {
    ## Marked point pattern
    marx <- marks(y, dfok=TRUE)
    if(is.data.frame(marx)) {
      ## select column or take first colum
      marx <- marx[, which.marks]
      y <- setmarks(y, marx)
    }
    if(npoints(y) > 0) {
      ok <- complete.cases(as.data.frame(y))
      if(!any(ok)) {
        warning("All mark values are NA; plotting locations only.")
        if(is.null(symap))
          symap <- default.symap.points(unmark(x), ..., chars=chars, cols=cols)
      } else if(any(!ok)) {
        warning(paste("Some marks are NA;",
                      "corresponding points are omitted."))
        x <- x[ok]
        y <- y[ok]
        marx <- marks(y)
      }
    }
    ## apply default symbol map
    if(is.null(symap))
      symap <- default.symap.points(y, chars=chars, cols=cols, 
                                    maxsize=maxsize, meansize=meansize,
                                    markscale=markscale,
                                    ...)
  }
#  gtype <- symbolmaptype(symap)

  ## Determine bounding box for main plot
  BB <- as.rectangle(x)
  sick <- inherits(x, "ppp") && !is.null(rejects <- attr(x, "rejects"))
  if(sick) {
    ## Get relevant parameters
    par.direct <- list(main=main, use.marks=use.marks,
                   maxsize=maxsize, meansize=meansize, markscale=markscale)
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
    sizeguess <- if(maxsize <= 0) NULL else (1.5 * maxsize)
    leg.args <- append(list(side=leg.side, vertical=vertical), leg.args)
    ## draw up layout
    legbox <- do.call.matched(plan.legend.layout,
                              append(list(B=BB, size = sizeguess,
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
  pt <- prepareTitle(main)
  main <- pt$main
  nlines <- pt$nlines
  blankmain <- if(nlines == 0) "" else rep("  ", nlines)
  cex.main <- resolve.1.default(list(cex.main=1), list(...))
  plot(BB, type="n", add=add, main=blankmain, show.all=show.all,
       cex.main=cex.main)

  if(sick) {
    if(show.window) {
      ## plot windows
      if(!is.null(rw)) {
        ## plot window for rejects
        rwinpardefault <- list(lty=2,lwd=1,border=1)
        rwinpars <-
          resolve.defaults(par.rejects, rwinpardefault)[names(rwinpardefault)]
        do.call(plot.owin, append(list(rwin, add=TRUE), rwinpars))
      }
      ## plot window of main pattern
      if(!clipped) {
        do.call(plot.owin,
                resolve.defaults(list(x$window, add=TRUE),
                                 list(...),
                                 list(invert=TRUE)))
      } else plot(clippy, add=TRUE, ...)
    }
    if(type != "n") {
      ## plot reject points
      do.call(plot.ppp, append(list(rejects, add=TRUE), par.all))
      warning(paste(rejects$n, "illegal points also plotted"))
    }
    ## the rest is added
    add <- TRUE
  }

  ## Now convert to bona fide point pattern
  x <- as.ppp(x)
  xwindow <- x$window

  ## Plot observation window (or at least the main title)
  do.call(plot.owin,
          resolve.defaults(list(x=xwindow,
                                add=TRUE,
                                main=main,
                                type=if(show.window && !clipped) "w" else "n",
                                show.all=show.all),
                           list(...),
                           list(invert=TRUE)))
  ## If clipped, plot visible part of original window
  if(show.window && clipped)
    plot(clippy, add=TRUE, ...)
  # else if(show.all) fakemaintitle(as.rectangle(xwindow), main, ...)

  if(type != "n") {
    ## plot symbols ##
    invoke.symbolmap(symap, marx, x, add=TRUE)
  }
  
  ## add legend
  if(legend) {
    b <- legbox$b
    legendmap <- if(length(leg.args) == 0) symap else 
                 do.call(update, append(list(object=symap), leg.args))
    do.call(plot,
            append(list(x=legendmap, main="", add=TRUE,
                        xlim=b$xrange, ylim=b$yrange),
                   leg.args))
  }
  
  return(invisible(symap))
}

plot.ppp

})


mark.scale.default <- function(marx, w, markscale=NULL,
                               maxsize=NULL, meansize=NULL,
                               characters=FALSE) {
  ## establish values of markscale, maxsize, meansize
  ngiven <- (!is.null(markscale)) +
            (!is.null(maxsize)) +
            (!is.null(meansize))
  if(ngiven > 1)
     stop("Only one of the arguments markscale, maxsize, meansize",
          " should be given", call.=FALSE)
  if(ngiven == 0) {
    ## if ALL are absent, enforce the spatstat defaults
    ## (which could also be null)
    pop <- spatstat.options("par.points")
    markscale <- pop$markscale
    maxsize   <- pop$maxsize
    meansize <- pop$meansize
  }
  ## Now check whether markscale is fixed
  if(!is.null(markscale)) {
    stopifnot(markscale > 0)
    return(markscale)
  }
  # Usual case: markscale is to be determined from maximum/mean physical size
  if(is.null(maxsize) && is.null(meansize)) {
    ## compute default value of 'maxsize'
    ## guess appropriate max physical size of symbols
    bb <- as.rectangle(w)
    maxsize <- 1.4/sqrt(pi * length(marx)/area(bb))
    maxsize <- min(maxsize, diameter(bb) * 0.07)
    ## updated: maxsize now represents *diameter*
    maxsize <- 2 * maxsize
  } else {
    if(!is.null(maxsize)) stopifnot(maxsize > 0) else stopifnot(meansize > 0)
  }
  
  # Examine mark values
  absmarx <- abs(marx)
  maxabs <- max(absmarx)
  tiny <- (maxabs < 4 * .Machine$double.eps)
  if(tiny)
    return(NA)

  ## finally determine physical scale for symbols
  if(!is.null(maxsize)) {
    scal <- maxsize/maxabs
  } else {
    meanabs <- mean(absmarx)
    scal <- meansize/meanabs
  }
  if(!characters) return(scal)

  ## if using characters ('pch') we need to
  ## convert physical sizes to 'cex' values
  charsize <- max(sidelengths(as.rectangle(w)))/40
  return(scal/charsize)
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
  do.call.matched(text.default,
                  resolve.defaults(list(x=x0, y=y0, labels=main),
                                   parlist,    list(...)),
                  funargs=graphicsPars("text"))
  return(invisible(NULL))
}
