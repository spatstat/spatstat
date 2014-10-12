#
#       plot.fv.R   (was: conspire.S)
#
#  $Revision: 1.116 $    $Date: 2014/10/12 08:17:57 $
#
#

conspire <- function(...) {
  .Deprecated("plot.fv", package="spatstat")
  plot.fv(...)
}

plot.fv <- local({

  hasonlyone <- function(x, amongst) {
    sum(all.vars(parse(text=x)) %in% amongst) == 1
  }

  extendifvector <- function(a, n, nmore) {
    if(is.null(a)) return(a)
    if(length(a) == 1) return(a)
    return(c(a, rep(a[1], nmore)))
  }

  fixit <- function(a, n, a0, a00) {
    # 'a' is formal argument
    # 'a0' and 'a00' are default and fallback default
    # 'n' is number of values required
    if(is.null(a))
      a <- if(!is.null(a0)) a0 else a00
    if(length(a) == 1)
      return(rep.int(a, n))
    else if(length(a) != n)
      stop(paste("Length of", short.deparse(substitute(a)),
                 "does not match number of curves to be plotted"))
    else 
      return(a)
  }

  pow10 <- function(x) { 10^x }

  plot.fv <- function(x, fmla, ..., subset=NULL, lty=NULL, col=NULL, lwd=NULL,
                      xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL,
                      ylim.covers=NULL, legend=!add, legendpos="topleft",
                      legendavoid=missing(legendpos),
                      legendmath=TRUE, legendargs=list(),
                      shade=fvnames(x, ".s"), shadecol="grey", add=FALSE,
                      log="",
                      mathfont=c("italic", "plain", "bold", "bolditalic"), 
                      limitsonly=FALSE) {

    xname <-
      if(is.language(substitute(x))) short.deparse(substitute(x)) else ""

    force(legendavoid)
    if(is.null(legend))
      legend <- !add

    mathfont <- match.arg(mathfont)

    verifyclass(x, "fv")
    env.user <- parent.frame()

    indata <- as.data.frame(x)

    xlogscale <- (log %in% c("x", "xy", "yx"))
    ylogscale <- (log %in% c("y", "xy", "yx"))

    ## ---------------- determine plot formula ----------------
  
    defaultplot <- missing(fmla) || is.null(fmla)
    if(defaultplot) 
      fmla <- formula(x)

    ## This *is* the last possible moment, so...
    fmla <- as.formula(fmla, env=env.user)

    ## validate the variable names
    vars <- variablesinformula(fmla)
    reserved <- c(".", ".x", ".y", ".a", ".s")
    external <- !(vars %in% c(colnames(x), reserved))
    if(any(external)) {
      sought <- vars[external]
      found <- unlist(lapply(sought, exists, envir=env.user, mode="numeric"))
      if(any(!found)) {
        nnot <- sum(!found)
        stop(paste(ngettext(nnot, "Variable", "Variables"),
                   commasep(sQuote(sought[!found])),
                   ngettext(nnot, "was", "were"),
                   "not found"))
      } else {
        ## validate the found variables
        externvars <- lapply(sought, get, envir=env.user)
        isnum <- unlist(lapply(externvars, is.numeric))
        len <- unlist(lapply(externvars, length))
        ok <- isnum & (len == 1 | len == nrow(x))
        if(!all(ok)) {
          nnot <- sum(!ok)
          stop(paste(ngettext(nnot, "Variable", "Variables"),
                     commasep(sQuote(sought[!ok])),
                     ngettext(nnot, "is", "are"),
                     "not of the right format"))
        }
      }
    }
  
    ## Extract left hand side as given
    lhs.original <- fmla[[2]]
    fmla.original <- fmla
  
    ## expand "."
    dotnames <- fvnames(x, ".")
    starnames <- fvnames(x, "*")
    umap <- fvexprmap(x)
    fmla <- eval(substitute(substitute(fom, um), list(fom=fmla, um=umap)))

    ## ------------------- extract data for plot ---------------------
  
    ## extract LHS and RHS of formula
    lhs <- fmla[[2]]
    rhs <- fmla[[3]]

    ## extract data 
    lhsdata <- eval(lhs, envir=indata)
    rhsdata <- eval(rhs, envir=indata)

    ## reformat
    if(is.vector(lhsdata)) {
      lhsdata <- matrix(lhsdata, ncol=1)
      lhsvars <- all.vars(as.expression(lhs))
      lhsvars <- lhsvars[lhsvars %in% names(x)]
      colnames(lhsdata) <-
        if(length(lhsvars) == 1) lhsvars else
        if(length(starnames) == 1 && (starnames %in% lhsvars)) starnames else 
        paste(short.deparse(lhs), collapse="")
    }
    ## check lhs names exist
    lnames <- colnames(lhsdata)
    nc <- ncol(lhsdata)
    lnames0 <- paste("V", seq_len(nc), sep="")
    if(length(lnames) != nc)
      colnames(lhsdata) <- lnames0
    else if(any(uhoh <- !nzchar(lnames)))
      colnames(lhsdata)[uhoh] <- lnames0[uhoh]
    lhs.names <- colnames(lhsdata)

    ## check whether each lhs column is associated with a single column of 'x'
    ## that is one of the alternative versions of the function.
    ##    This may be unreliable, as it depends on the
    ##    column names assigned to lhsdata by eval()
    one.star <- unlist(lapply(lhs.names, hasonlyone, amongst=fvnames(x, "*")))
    one.dot  <- unlist(lapply(lhs.names, hasonlyone, amongst=dotnames))
    explicit.lhs.names    <- ifelse(one.star, lhs.names, "")
    explicit.lhs.dotnames <- ifelse(one.star & one.dot, lhs.names, "")
  
    ## check rhs data
    if(is.matrix(rhsdata))
      stop("rhs of formula should yield a vector")
    rhsdata <- as.numeric(rhsdata)

    nplots <- ncol(lhsdata)
    allind <- 1:nplots
  
    ## ---------- extra plots may be implied by 'shade' -----------------
    extrashadevars <- NULL
  
    if(!is.null(shade)) {
      ## select columns by name or number
      names(allind) <- explicit.lhs.names
      shind <- try(allind[shade])
      if(inherits(shind, "try-error")) 
        stop(paste("The argument shade should be a valid subset index",
                   "for columns of x"), call.=FALSE)
      if(any(nbg <- is.na(shind))) {
        ## columns not included in formula: add them
        morelhs <- try(as.matrix(indata[ , shade[nbg], drop=FALSE]))
        if(inherits(morelhs, "try-error")) 
          stop(paste("The argument shade should be a valid subset index",
                     "for columns of x"), call.=FALSE)
        nmore <- ncol(morelhs)
        extrashadevars <- colnames(morelhs)
        if(defaultplot) {
          success <- TRUE
        } else if("." %in% variablesinformula(fmla.original)) {
          ## evaluate lhs of formula, expanding "." to shade names
          u <- if(length(extrashadevars) == 1) as.name(extrashadevars) else {
            as.call(lapply(c("cbind", extrashadevars), as.name))
          }
          ux <- as.name(fvnames(x, ".x"))
          uy <- as.name(fvnames(x, ".y"))
          foo <- eval(substitute(substitute(fom, list(.=u, .x=ux, .y=uy)),
                                 list(fom=fmla.original)))
          lhsnew <- foo[[2]]
          morelhs <- eval(lhsnew, envir=indata)
          success <- identical(colnames(morelhs), extrashadevars)
        } else if(is.name(lhs) && as.character(lhs) %in% names(indata)) {
          ## lhs is the name of a single column in x
          ## expand the LHS 
          explicit.lhs.names <- c(explicit.lhs.names, extrashadevars)
          ff <- paste("cbind",
                      paren(paste(explicit.lhs.names, collapse=", ")),
                      "~ 1")
          lhs <- lhs.of.formula(as.formula(ff))
          success <- TRUE
        } else if(length(explicit.lhs.dotnames) > 1) {
          ## lhs = cbind(...) where ... are dotnames
          cbound <- paste0("cbind",
                           paren(paste(explicit.lhs.dotnames, collapse=", ")))
          if(identical(deparse(lhs), cbound)) {
            success <- TRUE
            explicit.lhs.names <- union(explicit.lhs.names, extrashadevars)
            ff <- paste("cbind",
                        paren(paste(explicit.lhs.names, collapse=", ")),
                        "~ 1")
            lhs <- lhs.of.formula(as.formula(ff))
          } else success <- FALSE
        } else success <- FALSE
        if(success) {
          ## add these columns to the plotting data
          lhsdata <- cbind(lhsdata, morelhs)
          shind[nbg] <- nplots + seq_len(nmore)
          lty <- extendifvector(lty, nplots, nmore)
          col <- extendifvector(col, nplots, nmore)
          lwd <- extendifvector(lwd, nplots, nmore)
          nplots <- nplots + nmore
          ## update the names
          one.star <- unlist(lapply(explicit.lhs.names,
                                    hasonlyone, amongst=fvnames(x, "*")))
          one.dot  <- unlist(lapply(explicit.lhs.names,
                                    hasonlyone, amongst=dotnames))
          explicit.lhs.names    <- ifelse(one.star, explicit.lhs.names, "")
          explicit.lhs.dotnames <- ifelse(one.star & one.dot,
                                          explicit.lhs.names, "")
        } else {
          ## cannot add columns
          warning(paste("Shade",
                        ngettext(sum(nbg), "column", "columns"),
                        commasep(sQuote(shade[nbg])),
                        "were missing from the plot formula, and were omitted"))
          shade <- NULL
          extrashadevars <- NULL
        }
      }
    }

    ## -------------------- determine plotting limits ----------------------
  
    ## restrict data to subset if desired
    if(!is.null(subset)) {
      keep <- if(is.character(subset)) {
                eval(parse(text=subset), envir=indata)
              } else eval(subset, envir=indata)
      lhsdata <- lhsdata[keep, , drop=FALSE]
      rhsdata <- rhsdata[keep]
    }
  
    ## determine x and y limits and clip data to these limits
    if(is.null(xlim) && add) {
      ## x limits are determined by existing plot
      xlim <- par("usr")[1:2]
    }
    if(!is.null(xlim)) {
      ok <- !is.finite(rhsdata) | (xlim[1] <= rhsdata & rhsdata <= xlim[2])
      rhsdata <- rhsdata[ok]
      lhsdata <- lhsdata[ok, , drop=FALSE]
    } else {
      ## if we're using the default argument, use its recommended range
      if(rhs == fvnames(x, ".x")) {
        xlim <- attr(x, "alim") %orifnull% range(rhsdata, finite=TRUE)
        if(xlogscale && xlim[1] <= 0) 
          xlim[1] <- min(rhsdata[is.finite(rhsdata) & rhsdata > 0], na.rm=TRUE)
        ok <- !is.finite(rhsdata) | (rhsdata >= xlim[1] & rhsdata <= xlim[2])
        rhsdata <- rhsdata[ok]
        lhsdata <- lhsdata[ok, , drop=FALSE]
      } else { ## actual range of values to be plotted
        if(xlogscale) {
          ok <- is.finite(rhsdata) & (rhsdata > 0) & apply(lhsdata > 0, 1, any)
          xlim <- range(rhsdata[ok])
        } else {
          xlim <- range(rhsdata, na.rm=TRUE)
        }
      }
    }

    if(is.null(ylim)) {
      yok <- is.finite(lhsdata)
      if(ylogscale)
        yok <- yok & (lhsdata > 0)
      ylim <- range(lhsdata[yok],na.rm=TRUE)
    }
    if(!is.null(ylim.covers))
      ylim <- range(ylim, ylim.covers)

    ## return x, y limits only?
    if(limitsonly)
      return(list(xlim=xlim, ylim=ylim))
  
    ## -------------  work out how to label the plot --------------------

    ## extract plot labels, substituting function name
    labl <- fvlabels(x, expand=TRUE)
    ## create plot label map (key -> algebraic expression)
    map <- fvlabelmap(x) 

    ## ......... label for x axis ..................

    if(is.null(xlab)) {
      argname <- fvnames(x, ".x")
      if(as.character(fmla)[3] == argname) {
        ## The x axis variable is the default function argument.
        ArgString <- fvlabels(x, expand=TRUE)[[argname]]
        xexpr <- parse(text=ArgString)
        ## use specified font
        xexpr <- fontify(xexpr, mathfont)
        ## Add name of unit of length?
        ax <- summary(unitname(x))$axis
        if(is.null(ax)) {
          xlab <- xexpr
        } else {
          xlab <- expression(VAR ~ COMMENT)
          xlab[[1]][[2]] <- xexpr[[1]]
          xlab[[1]][[3]] <- ax
        }
      } else {
        ## map ident to label
        xlab <- eval(substitute(substitute(rh, mp), list(rh=rhs, mp=map)))
        ## use specified font
        xlab <- fontify(xlab, mathfont)
      }
    }
    if(is.language(xlab) && !is.expression(xlab))
      xlab <- as.expression(xlab)

    ## ......... label for y axis ...................

    leftside <- lhs
    if(ncol(lhsdata) > 1 || length(dotnames) == 1) {
      ## For labelling purposes only, simplify the LHS by 
      ## replacing 'cbind(.....)' by '.'
      ## even if not all columns are included.
      leftside <- paste(as.expression(leftside))
      eln <- explicit.lhs.dotnames
      eln <- eln[nzchar(eln)]
      cb <- if(length(eln) == 1) eln else {
        paste("cbind(",
              paste(eln, collapse=", "),
              ")", sep="")
      }
      compactleftside <- gsub(cb, ".", leftside, fixed=TRUE)
      ## Separately expand "." to cbind(.....)
      ## and ".x", ".y" to their real names
      dotdot <- c(dotnames, extrashadevars)
      cball <- if(length(dotdot) == 1) dotdot else {
        paste("cbind(",
              paste(dotdot, collapse=", "),
              ")", sep="")
      }
      expandleftside <- gsub(".x", fvnames(x, ".x"), leftside, fixed=TRUE)
      expandleftside <- gsub(".y", fvnames(x, ".y"), expandleftside, fixed=TRUE)
      expandleftside <- gsubdot(cball, expandleftside)
      ## convert back to language
      compactleftside <- parse(text=compactleftside)[[1]]
      expandleftside <- parse(text=expandleftside)[[1]]
    } else {
      compactleftside <- expandleftside <- leftside
    }

    ## construct label for y axis
    if(is.null(ylab)) {
      yl <- attr(x, "yexp")
      if(defaultplot && !is.null(yl)) {
        ylab <- yl
      } else {
        ## replace "." and short identifiers by plot labels
        ylab <- eval(substitute(substitute(le, mp),
                                list(le=compactleftside, mp=map)))
      }
    }
    if(is.language(ylab)) {
      ## use specified font
      ylab <- fontify(ylab, mathfont)
      ## ensure it's an expression
      if(!is.expression(ylab))
        ylab <- as.expression(ylab)
    }

    ## ------------------ start plotting ---------------------------

    ## create new plot
    if(!add)
      do.call("plot.default",
              resolve.defaults(list(xlim, ylim, type="n", log=log),
                               list(xlab=xlab, ylab=ylab),
                               list(...),
                               list(main=xname)))

    ## handle 'type' = "n" 
    giventype <- resolve.defaults(list(...), list(type=NA))$type
    if(identical(giventype, "n"))
      return(invisible(NULL))

    ## process lty, col, lwd arguments

    opt0 <- spatstat.options("par.fv")
  
    lty <- fixit(lty, nplots, opt0$lty, 1:nplots)
    col <- fixit(col, nplots, opt0$col, 1:nplots)
    lwd <- fixit(lwd, nplots, opt0$lwd, 1)

    ## convert to greyscale?
    if(spatstat.options("monochrome"))
      col <- to.grey(col)
    
    if(!is.null(shade)) {
      ## shade region between critical boundaries
      ## extract relevant columns for shaded bands
      shdata <- lhsdata[, shind]
      if(!is.matrix(shdata) || ncol(shdata) != 2) 
        stop("The argument shade should select two columns of x")
      ## truncate infinite values to plot limits
      if(any(isinf <- is.infinite(shdata))) {
        if(is.null(ylim)) {
          warning("Unable to truncate infinite values to the plot area")
        } else {
          shdata[isinf & (shdata == Inf)] <- ylim[2]
          shdata[isinf & (shdata == -Inf)] <- ylim[1]
        }
      }
      ## determine limits of shading
      shdata1 <- shdata[,1]
      shdata2 <- shdata[,2]
      ## plot grey polygon
      xpoly <- c(rhsdata, rev(rhsdata))
      ypoly <- c(shdata1, rev(shdata2)) 
      miss1 <- !is.finite(shdata1)
      miss2 <- !is.finite(shdata2)
      if(!any(broken <- (miss1 | miss2))) {
        ## single polygon
        polygon(xpoly, ypoly, border=shadecol, col=shadecol)
      } else {
        ## interrupted
        dat <- data.frame(rhsdata=rhsdata, shdata1=shdata1, shdata2=shdata2)
        serial <- cumsum(broken)
        lapply(split(dat, serial),
               function(z) {
                 with(z, {
                   xp <- c(rhsdata, rev(rhsdata))
                   yp <- c(shdata1, rev(shdata2))
                   polygon(xp, yp, border=shadecol, col=shadecol)
                 })
               })
        ## save for use in placing legend
        okp <- !c(broken, rev(broken))
        xpoly <- xpoly[okp]
        ypoly <- ypoly[okp]
      }
      ## overwrite graphical parameters
      lty[shind] <- 1
      ## try to preserve the same type of colour specification
      if(is.character(col) && is.character(shadecol)) {
        ## character representations 
        col[shind] <- shadecol
      } else if(is.numeric(col) && !is.na(sc <- paletteindex(shadecol))) {
        ## indices in colour palette
        col[shind] <- sc
      } else {
        ## convert colours to hexadecimal and edit relevant values
        col <- col2hex(col)
        col[shind] <- col2hex(shadecol)
      }
      ## remove these columns from further plotting
      allind <- allind[-shind]
      ## 
    } else xpoly <- ypoly <- numeric(0)
  
    ## ----------------- plot lines ------------------------------

    for(i in allind)
      lines(rhsdata, lhsdata[,i], lty=lty[i], col=col[i], lwd=lwd[i])

    if(nplots == 1)
      return(invisible(NULL))

    ## ---------------- determine legend -------------------------
    key <- colnames(lhsdata)
    mat <- match(key, names(x))
    keyok <- !is.na(mat)
    matok <- mat[keyok]
    legdesc <- rep.int("constructed variable", length(key))
    legdesc[keyok] <- attr(x, "desc")[matok]
    leglabl <- lnames0
    leglabl[keyok] <- labl[matok]
    ylab <- attr(x, "ylab")
    if(!is.null(ylab)) {
      if(is.language(ylab)) 
        ylab <- flat.deparse(ylab)
      legdesc <- sprintf(legdesc, ylab)
    }
    ## compute legend info
    legtxt <- key
    if(legendmath) {
      legtxt <- leglabl
      if(defaultplot) {
        ## try to convert individual labels to expressions
        fancy <- try(parse(text=leglabl), silent=TRUE)
        if(!inherits(fancy, "try-error"))
          legtxt <- fancy
      } else {
        ## try to navigate the parse tree
        fancy <- try(fvlegend(x, expandleftside), silent=TRUE)
        if(!inherits(fancy, "try-error"))
          legtxt <- fancy
      }
    }

    if(is.expression(legtxt) ||
       is.language(legtxt) ||
       all(sapply(legtxt, is.language)))
      legtxt <- fontify(legtxt, mathfont)

    ## --------------- handle legend plotting  -----------------------------
    
    if(identical(legend, TRUE)) {
      ## legend will be plotted
      ## Basic parameters of legend
      legendxpref <- if(identical(legendpos, "float")) NULL else legendpos
      optparfv <- spatstat.options("par.fv")$legendargs %orifnull% list()
      legendspec <- resolve.defaults(legendargs,
                                     list(lty=lty,
                                          col=col,
                                          lwd=lwd),
                                     optparfv,
                                     list(x=legendxpref,
                                          legend=legtxt,
                                          inset=0.05,
                                          y.intersp=if(legendmath) 1.3 else 1),
                                     .StripNull=TRUE)
      
      if(legendavoid || identical(legendpos, "float")) {
        ## Automatic determination of legend position
        ## Assemble data for all plot objects
        linedata <- list()
        xmap <- if(xlogscale) log10 else identity
        ymap <- if(ylogscale) log10 else identity
        inv.xmap <- if(xlogscale) pow10 else identity
        inv.ymap <- if(ylogscale) pow10 else identity 
        for(i in seq_along(allind)) 
          linedata[[i]] <- list(x=xmap(rhsdata), y=ymap(lhsdata[,i]))
        polydata <-
          if(length(xpoly) > 0) list(x=xmap(xpoly), y=ymap(ypoly)) else NULL
        objects <- assemble.plot.objects(xmap(xlim), ymap(ylim),
                                         lines=linedata, polygon=polydata)
        ## find best position to avoid them
        legendbest <- findbestlegendpos(objects, preference=legendpos,
                                      legendspec=legendspec)
        ## handle log scale
        if((xlogscale || ylogscale) &&
           checkfields(legendbest, c("x", "xjust", "yjust"))) {
          ## back-transform x, y coordinates
          legendbest$x$x <- inv.xmap(legendbest$x$x)
          legendbest$x$y <- inv.ymap(legendbest$x$y)
        }
      } else legendbest <- list()
    
      ##  ********** plot legend *************************
      if(!is.null(legend) && legend) 
        do.call("legend",
                resolve.defaults(legendargs,
                                 legendbest,
                                 legendspec,
                                 .StripNull=TRUE))
      
    }

    ## convert labels back to character
    labl <- paste.expr(legtxt)
    labl <- gsub(" ", "", labl)
    ## return legend info
    df <- data.frame(lty=lty, col=col, key=key, label=labl,
                     meaning=legdesc, row.names=key)
    return(df)
  }
  plot.fv

})



assemble.plot.objects <- function(xlim, ylim, ..., lines=NULL, polygon=NULL) {
  # Take data that would have been passed to the commands 'lines' and 'polygon'
  # and form corresponding geometrical objects.
  objects <- list()
  if(!is.null(lines)) {
    if(is.psp(lines)) {
      objects <- list(lines)
    } else {
      if(checkfields(lines, c("x", "y"))) {
        lines <- list(lines)
      } else if(!all(unlist(lapply(lines, checkfields, L=c("x", "y")))))
        stop("lines should be a psp object, a list(x,y) or a list of list(x,y)")
      W <- owin(xlim, ylim)
      for(i in seq_along(lines)) {
        lines.i <- lines[[i]]
        x.i <- lines.i$x
        y.i <- lines.i$y
        n <- length(x.i)
        if(length(y.i) != n)
          stop(paste(paste("In lines[[", i, "]]", sep=""),
                     "the vectors x and y have unequal length"))
        if(!all(ok <- (is.finite(x.i) & is.finite(y.i)))) {
          x.i <- x.i[ok]
          y.i <- y.i[ok]
          n <- sum(ok)
        }
        segs.i <- psp(x.i[-n], y.i[-n], x.i[-1], y.i[-1], W, check=FALSE)
        objects <- append(objects, list(segs.i))        
      }
    }
  }
  if(!is.null(polygon)) {
    # Add filled polygon
    pol <- polygon[c("x", "y")]
    ok <- with(pol, is.finite(x) & is.finite(y))
    if(!all(ok))
      pol <- with(pol, list(x=x[ok], y=y[ok]))
    if(Area.xypolygon(pol) < 0) pol <- lapply(pol, rev)
    P <- try(owin(poly=pol, xrange=xlim, yrange=ylim, check=FALSE))
    if(!inherits(P, "try-error"))
      objects <- append(objects, list(P))
  }
  return(objects)
}

findbestlegendpos <- local({
  # Given a list of geometrical objects, find the best position
  # to avoid them.
  thefunction <- function(objects, show=FALSE, aspect=1, bdryok=TRUE,
                          preference="float", verbose=FALSE,
                          legendspec=NULL) {
    # find bounding box
    W <- do.call("boundingbox", lapply(objects, as.rectangle))
    # convert to common box
    objects <- lapply(objects, rebound, rect=W)
    # comp
    # rescale x and y axes so that bounding box has aspect ratio 'aspect'
    aspectW <- with(W, diff(yrange)/diff(xrange))
    s <- aspect/aspectW
    mat <- diag(c(1, s))
    invmat <- diag(c(1, 1/s))
    scaled.objects <- lapply(objects, affine, mat=mat)
    scaledW <- affine(W, mat=mat)
    if(verbose) {
      cat("Scaled space:\n")
      print(scaledW)
    }
    # pixellate the scaled objects
    pix.scal.objects <- lapply(scaled.objects, asma)
    # apply distance transforms in scaled space
    D1 <- distmap(pix.scal.objects[[1]])
    Dlist <- lapply(pix.scal.objects, distmap, xy=list(x=D1$xcol, y=D1$yrow))
    # distance transform of superposition
    D <- im.apply(Dlist, min)
    if(!bdryok) {
      # include distance to boundary
      B <- attr(D1, "bdry")
      D <- eval.im(pmin.int(D, B))
    }
    if(show) {
      plot(affine(D, mat=invmat), add=TRUE)
      lapply(lapply(scaled.objects, affine, mat=invmat), plot, add=TRUE)
    }
    if(preference != "float") {
      # evaluate preferred location (check for collision)
      if(!is.null(legendspec)) {
        # pretend to plot the legend as specified
        legout <- do.call("legend", append(legendspec, list(plot=FALSE)))
        # determine bounding box
        legbox <- with(legout$rect, owin(c(left, left+w), c(top-h, top)))
        scaledlegbox <- affine(legbox, mat=mat)
        # check for collision 
        Dmin <- min(D[scaledlegbox])
        if(Dmin >= 0.02) {
          # no collision: stay at preferred location. Exit.
          return(list(x=preference))
        }
        # collision occurred! 
      } else {
        # no legend information.
        # Pretend legend is 15% of plot width and height
        xr <- scaledW$xrange
        yr <- scaledW$yrange
        testloc <- switch(preference,
                          topleft     = c(xr[1],yr[2]),
                          top         = c(mean(xr), yr[2]),
                          topright    = c(xr[2], yr[2]),
                          right       = c(xr[2], mean(yr)),
                          bottomright = c(xr[2], yr[1]),
                          bottom      = c(mean(xr), yr[1]),
                          bottomleft  = c(xr[1], yr[1]),
                          left        = c(xr[1], mean(yr)),
                          center      = c(mean(xr), mean(yr)),
                          NULL)
        if(!is.null(testloc)) {
          # look up distance value at preferred location
          val <- safelookup(D, list(x=testloc[1], y=testloc[2]))
          crit <- 0.15 * min(diff(xr), diff(yr))
          if(verbose)
            cat(paste("val=",val, ", crit=", crit, "\n"))
          if(val > crit) {
            # no collision: stay at preferred location. Exit.
            return(list(x=preference))
          }
        # collision occurred! 
        }
      }
      # collision occurred! 
    }
    # find location of max
    locmax <- which(D$v == max(D), arr.ind=TRUE)
    locmax <- unname(locmax[1,])
    pos <- list(x=D$xcol[locmax[2]], y=D$yrow[locmax[1]])
    pos <- affinexy(pos, mat=invmat)
    if(show) 
      points(pos)
    # determine justification of legend relative to this point
    # to avoid crossing edges of plot
    xrel <- (pos$x - W$xrange[1])/diff(W$xrange)
    yrel <- (pos$y - W$yrange[1])/diff(W$yrange)
    xjust <- if(xrel < 0.1) 0 else if(xrel > 0.9) 1 else 0.5 
    yjust <- if(yrel < 0.1) 0 else if(yrel > 0.9) 1 else 0.5
    #
    out <- list(x=pos, xjust=xjust, yjust=yjust)
    return(out)
  }

  asma <- function(z) { if(is.owin(z)) as.mask(z) else
                        if(is.psp(z)) as.mask.psp(z) else NULL }
  
  callit <- function(...) {
    rslt <- try(thefunction(...))
    if(!inherits(rslt, "try-error"))
      return(rslt)
    return(list())
  }
  callit
})
  
