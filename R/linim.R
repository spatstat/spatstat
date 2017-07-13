#
# linim.R
#
#  $Revision: 1.35 $   $Date: 2017/07/13 02:43:30 $
#
#  Image/function on a linear network
#

linim <- function(L, Z, ..., df=NULL) {
  L <- as.linnet(L)
  stopifnot(is.im(Z))
  if(is.null(df)) {
    # compute the data frame of mapping information
    xx <- rasterx.im(Z)
    yy <- rastery.im(Z)
    mm <- !is.na(Z$v)
    xx <- as.vector(xx[mm])
    yy <- as.vector(yy[mm])
    pixelcentres <- ppp(xx, yy, window=as.rectangle(Z), check=FALSE)
    pixdf <- data.frame(xc=xx, yc=yy)
    # project pixel centres onto lines
    p2s <- project2segment(pixelcentres, as.psp(L))
    projloc <- as.data.frame(p2s$Xproj)
    projmap <- as.data.frame(p2s[c("mapXY", "tp")])
    # extract values
    values <- Z[pixelcentres]
    # bundle
    df <- cbind(pixdf, projloc, projmap, data.frame(values=values))
  } else {
    stopifnot(is.data.frame(df))
    neednames <- c("xc", "yc", "x", "y", "mapXY", "tp", "values")
    ok <- neednames %in% names(df)
    if(any(!ok)) {
      nn <- sum(!ok)
      stop(paste(ngettext(nn, "A column", "Columns"),
                 "named", commasep(sQuote(neednames[!ok])),
                 ngettext(nn, "is", "are"),
                 "missing from argument", sQuote("df")))
    }
  }
  out <- Z
  attr(out, "L") <- L
  attr(out, "df") <- df
  class(out) <- c("linim", class(out))
  return(out)
}

print.linim <- function(x, ...) {
  splat("Image on linear network")
  print(attr(x, "L"))
  NextMethod("print")
}

summary.linim <- function(object, ...) {
  y <- NextMethod("summary")
  if("integral" %in% names(y))
    y$integral <- integral(object)
  y$network <- summary(as.linnet(object))
  class(y) <- c("summary.linim", class(y))
  return(y)
}

print.summary.linim <- function(x, ...) {
  splat(paste0(x$type, "-valued"), "pixel image on a linear network")
  unitinfo <- summary(x$units)
  pluralunits <- unitinfo$plural
  sigdig <- getOption('digits')
  di <- x$dim
  win <- x$window
  splat(di[1L], "x", di[2L], "pixel array (ny, nx)")
  splat("enclosing rectangle:",
        prange(signif(win$xrange, sigdig)),
        "x",
        prange(signif(win$yrange, sigdig)),
        unitinfo$plural,
        unitinfo$explain)
  splat("dimensions of each pixel:",
        signif(x$xstep, 3), "x", signif(x$ystep, sigdig),
        pluralunits)
  if(!is.null(explain <- unitinfo$explain))
    splat(explain)
  splat("Pixel values (on network):")
  switch(x$type,
         integer=,
         real={
           splat("\trange =", prange(signif(x$range, sigdig)))
           splat("\tintegral =", signif(x$integral, sigdig))
           splat("\tmean =", signif(x$mean, sigdig))
         },
         factor={
           print(x$table)
         },
         complex={
           splat("\trange: Real",
                 prange(signif(x$Re$range, sigdig)),
                 "Imaginary",
                 prange(signif(x$Im$range, sigdig)))
           splat("\tintegral =", signif(x$integral, sigdig))
           splat("\tmean =", signif(x$mean, sigdig))
         },
         {
           print(x$summary)
         })
  splat("Underlying network:")
  print(x$network)
  return(invisible(NULL))
}


plot.linim <- function(x, ..., style=c("colour", "width"),
                       scale, adjust=1,
		       legend=TRUE,
                       leg.side=c("right", "left", "bottom", "top"),
                       leg.sep=0.1,
                       leg.wid=0.1,
                       leg.args=list(),
                       leg.scale=1,
                       do.plot=TRUE) {
  xname <- short.deparse(substitute(x))
  style <- match.arg(style)
  leg.side <- match.arg(leg.side)
  ribstuff <- list(ribside  = leg.side,
                   ribsep   = leg.sep,
                   ribwid   = leg.wid,
                   ribargs  = leg.args,
                   ribscale = leg.scale)
  # colour style: plot as pixel image
  if(style == "colour" || !do.plot)
    return(do.call(plot.im,
                   resolve.defaults(list(x),
                                    list(...),
                                    ribstuff,
                                    list(main=xname,
				         legend=legend,
					 do.plot=do.plot))))
  # width style
  L <- attr(x, "L")
  df <- attr(x, "df")
  Llines <- as.psp(L)
  W <- as.owin(L)
  # plan layout
  if(legend) {
    # use layout procedure in plot.im
    z <- do.call(plot.im,
		 resolve.defaults(list(x, do.plot=FALSE, legend=TRUE),
                                  list(...),
                                  ribstuff,
                                  list(main=xname)))
    bb.all <- attr(z, "bbox")
    bb.leg <- attr(z, "bbox.legend")
  } else {
    bb.all <- Frame(W)
    bb.leg <- NULL
  }
  legend <- !is.null(bb.leg)
  if(legend) {
    # expand plot region to accommodate text annotation in legend
    if(leg.side %in% c("left", "right")) {
      delta <- 2 * sidelengths(bb.leg)[1]
      xmargin <- if(leg.side == "right") c(0, delta) else c(delta, 0)
      bb.all <- grow.rectangle(bb.all, xmargin=xmargin)
    }
  }
  # initialise plot
  bb <- do.call.matched(plot.owin,
                        resolve.defaults(list(x=bb.all, type="n"),
                                         list(...), list(main=xname)),
                        extrargs="type")
  # resolve graphics parameters for polygons
  grafpar <- resolve.defaults(list(...), list(border=1, col=1))
  grafpar <- grafpar[names(grafpar) %in% names(formals(polygon))]
  # rescale values to a plottable range
  vr <- range(df$values)
  vr[1L] <- min(0, vr[1L])
  if(missing(scale)) {
    maxsize <- mean(distmap(Llines))/2
    scale <- maxsize/diff(vr)
  } 
  df$values <- adjust * scale * (df$values - vr[1L])/2
  # split data by segment
  mapXY <- factor(df$mapXY, levels=seq_len(Llines$n))
  dfmap <- split(df, mapXY, drop=TRUE)
  # sort each segment's data by position along segment
  dfmap <- lapply(dfmap, sortalongsegment)
  # plot each segment's data
#  Lends <- Llines$ends
  Lperp <- angles.psp(Llines) + pi/2
  Lfrom <- L$from
  Lto   <- L$to
  Lvert <- L$vertices
  Ljoined  <- (vertexdegree(L) > 1)
  # precompute coordinates of dodecagon
  dodo <- disc(npoly=12)$bdry[[1L]]
  #
  for(i in seq(length(dfmap))) {
    z <- dfmap[[i]]
    segid <- unique(z$mapXY)[1L]
    xx <- z$x
    yy <- z$y
    vv <- z$values
    # add endpoints of segment
    ileft <- Lfrom[segid]
    iright <- Lto[segid]
    leftend <- Lvert[ileft]
    rightend <- Lvert[iright]
    xx <- c(leftend$x, xx, rightend$x)
    yy <- c(leftend$y, yy, rightend$y)
    vv <- c(vv[1L],     vv, vv[length(vv)])
    rleft <- vv[1L]
    rright <- vv[length(vv)]
    # draw polygon around segment
    xx <- c(xx, rev(xx))
    yy <- c(yy, rev(yy))
    vv <- c(vv, -rev(vv))
    ang <- Lperp[segid]
    xx <- xx + cos(ang) * vv
    yy <- yy + sin(ang) * vv
    ## first add dodecagonal 'joints'
    if(Ljoined[ileft] && rleft > 0) 
      do.call(polygon,
              append(list(x=rleft * dodo$x + leftend$x,
                          y=rleft * dodo$y + leftend$y),
                     grafpar))
    if(Ljoined[iright] && rright > 0)
      do.call(polygon,
              append(list(x=rright * dodo$x + rightend$x,
                          y=rright * dodo$y + rightend$y),
                     grafpar))
    # now draw main
    do.call(polygon, append(list(x=xx, y=yy), grafpar))
  }
  result <- adjust * scale
  attr(result, "bbox") <- bb
  if(legend) {
    attr(result, "bbox.legend") <- bb.leg
    ## get graphical arguments
    grafpar <- resolve.defaults(leg.args, grafpar)
    grafpar <- grafpar[names(grafpar) %in% names(formals(polygon))]
    ## set up scale of typical pixel values
    gvals <- leg.args$at %orifnull% prettyinside(range(x))
    # corresponding widths
    wvals <- adjust * scale * gvals
    # glyph positions
    ng <- length(gvals)
    xr <- bb.leg$xrange
    yr <- bb.leg$yrange
    switch(leg.side,
           right = ,
	   left = {
	     y <- seq(yr[1], yr[2], length.out=ng+1L)
	     y <- (y[-1L] + y[-(ng+1L)])/2
	     for(j in 1:ng) {
               xx <- xr[c(1L,2L,2L,1L)]
	       yy <- (y[j] + c(-1,1) * wvals[j]/2)[c(1L,1L,2L,2L)]
	       do.call(polygon, append(list(xx, yy), grafpar))
	     }
	   },
	   bottom = ,
	   top = {
	     x <- seq(xr[1], xr[2], length.out=ng+1L)
	     x <- (x[-1L] + x[-(ng+1L)])/2
	     for(j in 1:ng) {
	       xx <- (x[j] + c(-1,1) * wvals[j]/2)[c(1L,1L,2L,2L)]
               yy <- yr[c(1L,2L,2L,1L)]
	       do.call(polygon, append(list(xx, yy), grafpar))
	     }
	   })
     # add text labels
     check.1.real(leg.scale)
     glabs <- leg.args$labels %orifnull% signif(leg.scale * gvals, 2)
     switch(leg.side,
            right  = text(xr[2], y,     pos=4, labels=glabs),
            left   = text(xr[1], y,     pos=2, labels=glabs),
	    bottom = text(x,     yr[1], pos=1, labels=glabs),
	    top    = text(x,     yr[2], pos=3, labels=glabs))
  }
  return(invisible(result))
}

sortalongsegment <- function(df) {
  df[fave.order(df$tp), , drop=FALSE]
}

as.im.linim <- function(X, ...) {
  attr(X, "L") <- attr(X, "df") <- NULL
  class(X) <- "im"
  if(length(list(...)) > 0)
    X <- as.im(X, ...)
  return(X)
}

as.linim <- function(X, ...) {
  UseMethod("as.linim")
}

as.linim.default <- function(X, L, ..., eps = NULL, dimyx = NULL, xy = NULL,
                                        delta = NULL) {
  stopifnot(inherits(L, "linnet"))
  Y <- as.im(X, W=as.rectangle(as.owin(L)), ..., eps=eps, dimyx=dimyx, xy=xy)
  M <- as.mask.psp(as.psp(L), as.owin(Y))
  Y[complement.owin(M)] <- NA
  df <- NULL
  if(!is.null(delta)) {
    df <- pointsAlongNetwork(L, delta)
    pix <- nearest.valid.pixel(df$x, df$y, Y)
    df$xc <- Y$xcol[pix$col]
    df$yc <- Y$yrow[pix$row]
    df$values <- Y$v[cbind(pix$row, pix$col)]
    df <- df[,c("xc", "yc", "x", "y", "seg", "tp", "values")]
    names(df)[names(df) == "seg"] <- "mapXY"
  }
  out <- linim(L, Y, df=df)
  return(out)
}

pointsAlongNetwork <- local({

  pointsAlongNetwork <- function(L, delta) {
    #' sample points evenly spaced along each segment
    stopifnot(inherits(L, "linnet"))
    S <- as.psp(L)
    ns <- nsegments(S)
    seglen <- lengths.psp(S)
    ends <- as.data.frame(S)
    nsample <- pmax(1, ceiling(seglen/delta))
    df <- NULL
    x0 <- ends$x0
    y0 <- ends$y0
    x1 <- ends$x1
    y1 <- ends$y1
    for(i in seq_len(ns)) {
      nn <- nsample[i] + 1L
      tcut <- seq(0, 1, length.out=nn)
      tp <- (tcut[-1] + tcut[-nn])/2
      x <- x0[i] * (1-tp) + x1[i] * tp
      y <- y0[i] * (1-tp) + y1[i] * tp
      df <- rbind(df, data.frame(x=x, y=y, seg=i, tp=tp))
    }
    return(df)          
  }

  pointsAlongNetwork
})

as.linim.linim <- function(X, ...) {
  if(length(list(...)) == 0)
    return(X)
  Y <- as.linim.default(X, as.linnet(X), ...)
  return(Y)
}

# analogue of eval.im

eval.linim <- function(expr, envir, harmonize=TRUE) {
  sc <- sys.call()
  # Get names of all variables in the expression
  e <- as.expression(substitute(expr))
  varnames <- all.vars(e)
  allnames <- all.names(e, unique=TRUE)
  funnames <- allnames[!(allnames %in% varnames)]
  if(length(varnames) == 0)
    stop("No variables in this expression")
  # get the values of the variables
  if(missing(envir)) {
    envir <- parent.frame() # WAS: sys.parent()
  } else if(is.list(envir)) {
    envir <- list2env(envir, parent=parent.frame())
  }
  vars <- mget(varnames, envir=envir, inherits=TRUE, ifnotfound=list(NULL))
  funs <- mget(funnames, envir=envir, inherits=TRUE, ifnotfound=list(NULL))
  # Find out which variables are (linear) images
  islinim <- unlist(lapply(vars, inherits, what="linim"))
  if(!any(islinim))
    stop("There are no linear images (class linim) in this expression")
  # ....................................
  # Evaluate the pixel values using eval.im
  # ....................................
  sc[[1L]] <- as.name('eval.im')
  sc$envir <- envir
  Y <- eval(sc)
  # .........................................
  # Then evaluate data frame entries if feasible
  # .........................................
  dfY <- NULL
  linims <- vars[islinim]
  nlinims <- length(linims)
  dframes <- lapply(linims, attr, which="df")
  nets <- lapply(linims, attr, which="L")
  isim <- unlist(lapply(vars, is.im))
  if(!any(isim & !islinim)) {
    # all images are 'linim' objects
    # Check that the images refer to the same linear network
    if(nlinims > 1) {
      agree <- unlist(lapply(nets[-1L], identical, y=nets[[1L]]))
      if(!all(agree))
        stop(paste("Images do not refer to the same linear network"))
    }
    dfempty <- unlist(lapply(dframes, is.null))
    if(!any(dfempty)) {
      # ensure data frames are compatible
      if(length(dframes) > 1 && (
          length(unique(nr <- sapply(dframes, nrow))) > 1   ||
           !allElementsIdentical(dframes, "seg")            ||
   	   !allElementsIdentical(dframes, "tp")
	)) {
        # find the one with finest spacing
	imax <- which.max(nr)
	# resample the others
	dframes[-imax] <- lapply(dframes[-imax],
	                         resampleNetworkDataFrame,
	                         template=dframes[[imax]])
      }
      # replace each image variable by its data frame column of values
      vars[islinim] <- lapply(dframes, getElement, "values")
      # now evaluate expression
      Yvalues <- eval(e, append(vars, funs))
      # pack up
      dfY <- dframes[[1L]]
      dfY$values <- Yvalues
    }
  }
  result <- linim(nets[[1L]], Y, df=dfY)
  return(result)
}

resampleNetworkDataFrame <- function(df, template) {
  # resample 'df' at the points of 'template'
  invalues  <- df$values
  insegment <- df$mapXY
  inteepee  <- df$tp
  out <- template
  n <- nrow(out)
  outvalues <- vector(mode = typeof(invalues), length=n)
  outsegment <- out$mapXY
  outteepee  <- out$tp
  for(i in seq_len(n)) {
    relevant <- which(insegment == outsegment[i])
    if(length(relevant) > 0) {
      j <- which.min(abs(inteepee[relevant] - outteepee[i]))
      outvalues[i] <- invalues[relevant[j]]
    }
  }
  out$values <- outvalues
  return(out)
}

as.linnet.linim <- function(X, ...) {
  attr(X, "L")
}

"[.linim" <- function(x, i, ..., drop=TRUE) {
  if(!missing(i) && is.lpp(i)) {
    n <- npoints(i)
    result <- vector(mode=typeof(x$v), length=n)
    if(n == 0) return(result)
    if(!is.null(df <- attr(x, "df"))) {
      #' use data frame of sample points along network
      knownseg <- df$mapXY
      knowntp  <- df$tp
      knownval <- df$values
      #' extract local coordinates of query points
      coo <- coords(i)
      queryseg <- coo$seg
      querytp  <- coo$tp
      #' match to nearest sample point
      for(j in 1:n) {
        relevant <- (knownseg == queryseg[j])
        if(!any(relevant)) {
          result[j] <- NA
        } else {
          k <- which.min(abs(knowntp[relevant] - querytp[j]))
          result[j] <- knownval[relevant][k]
        }
      }
      if(drop && anyNA(result))
        result <- result(!is.na(result))
      return(result)
    }
    #' give up and use pixel image
  }
  #' apply subset method for 'im'
  y <- NextMethod("[")
  if(!is.im(y)) return(y) # vector of pixel values
  #' handle linear network info
  L <- attr(x, "L")
  df <- attr(x, "df")
  #' clip to new window
  W <- Window(y)
  attr(y, "L") <- L[W]
  attr(y, "df") <- df[inside.owin(df$xc, df$yc, W), , drop=FALSE]
  return(y)
}

integral.linim <- function(f, domain=NULL, ...){
  verifyclass(f, "linim")
  if(!is.null(domain)) 
    f <- f[domain]
  #' extract data
  L <- as.linnet(f)
  ns <- nsegments(L)
  df <- attr(f, "df")
  vals <- df$values
  seg <- factor(df$mapXY, levels=1:ns)
  #' ensure each segment has at least one sample point
  nper <- table(seg)
  if(any(missed <- (nper == 0))) {
    missed <- unname(which(missed))
    xy <- midpoints.psp(as.psp(L)[missed])
    valxy <- f[xy]
    seg <- c(seg, factor(missed, levels=1:ns))
    vals <- c(vals, valxy)
  }
  #' take average of data on each segment
  mu <- as.numeric(by(vals, seg, mean, ..., na.rm=TRUE))
  mu[is.na(mu)] <- 0
  #' weighted sum
  len <- lengths.psp(as.psp(L))
  if(anyNA(vals)) {
    p <- as.numeric(by(!is.na(vals), seg, mean, ..., na.rm=TRUE))
    p[is.na(p)] <- 0
    len <- len * p
  }
  return(sum(mu * len))
}

mean.linim <- function(x, ...) {
  trap.extra.arguments(...)
  integral(x)/sum(lengths.psp(as.psp(as.linnet(x))))
}

quantile.linim <- function(x, probs = seq(0,1,0.25), ...) {
  verifyclass(x, "linim")
  #' extract data
  df <- attr(x, "df")
  L <- as.linnet(x)
  vals <- df$values
  #' count sample points on each segment
  seg <- factor(df$mapXY, levels=1:nsegments(L))
  nvals <- table(seg)
  #' calculate weights
  len <- lengths.psp(as.psp(L))
  iseg <- as.integer(seg)
  wts <- len[iseg]/nvals[iseg]
  return(weighted.quantile(vals, wts, probs))
}

median.linim <- function(x, ...) {
  trap.extra.arguments(...)
  return(unname(quantile(x, 0.5)))
}

shift.linim <- function (X, ...) {
  verifyclass(X, "linim")
  Z <- shift(as.im(X), ...)
  L <- shift(as.linnet(X), ...)
  v <- getlastshift(L)
  df <- attr(X, "df")
  df[,c("xc","yc")] <- shiftxy(df[,c("xc", "yc")], v)
  df[,c("x","y")]   <- shiftxy(df[,c("x", "y")],   v)
  Y <- linim(L, Z, df=df)
  return(putlastshift(Y, v))
}

affine.linim <- function(X, mat = diag(c(1, 1)), vec = c(0, 0), ...) {
  Z <- affine(as.im(X), mat=mat, vec=vec, ...)
  L <- affine(as.linnet(X), mat=mat, vec=vec, ...)
  df <- attr(X, "df")
  df[,c("xc","yc")] <- affinexy(df[,c("xc", "yc")], mat=mat, vec=vec)
  df[,c("x","y")]   <- affinexy(df[,c("x", "y")],   mat=mat, vec=vec)
  Y <- linim(L, Z, df=df)
  return(Y)
}

scalardilate.linim <- function(X, f, ..., origin=NULL) {
  trap.extra.arguments(..., .Context = "In scalardilate(X,f)")
  check.1.real(f, "In scalardilate(X,f)")
  stopifnot(is.finite(f) && f > 0)
  if (!is.null(origin)) {
    X <- shift(X, origin = origin)
    negorig <- getlastshift(X)
  }
  else negorig <- c(0, 0)
  Y <- affine(X, mat = diag(c(f, f)), vec = -negorig)
  return(Y)
}

as.data.frame.linim <- function(x, ...) {
  df <- attr(x, "df")
  if(!is.na(m <- match("mapXY", colnames(df))))
    colnames(df)[m] <- "seg"
  return(df)
}

pairs.linim <- function(..., plot=TRUE, eps=NULL) {
  argh <- list(...)
  ## unpack single argument which is a list of images
  if(length(argh) == 1) {
    arg1 <- argh[[1L]]
    if(is.list(arg1) && all(sapply(arg1, is.im)))
      argh <- arg1
  }
  ## identify which arguments are images
  isim <- sapply(argh, is.im)
  nim <- sum(isim)
  if(nim == 0) 
    stop("No images provided")
  ## separate image arguments from others
  imlist <- argh[isim]
  rest   <- argh[!isim]
  ## identify which arguments are images on a network
  islinim <- sapply(imlist, inherits, what="linim")
  if(!any(islinim)) # shouldn't be here
    return(pairs.im(argh, plot=plot))
  ## adjust names
  imnames <- names(imlist) %orifnull% rep("", length(imlist))
  if(any(needsname <- !nzchar(imnames))) 
    imnames[needsname] <- paste0("V", seq_len(nim)[needsname])
  names(imlist) <- imnames
  ## choose resolution
  if(is.null(eps)) {
    xstep <- min(sapply(imlist, getElement, name="xstep"))
    ystep <- min(sapply(imlist, getElement, name="ystep"))
    eps <- min(xstep, ystep)
  }
  ## extract linear network
  Z1 <- imlist[[min(which(islinim))]]
  L <- as.linnet(Z1)
  ## construct equally-spaced sample points
  X <- pointsOnLines(as.psp(L), eps=eps)
  ## sample each image
  pixvals <- lapply(imlist, "[", i=X, drop=FALSE)
  pixdf <- as.data.frame(pixvals)
  ## pairs plot
  if(plot) {
    if(nim > 1) {
      do.call(pairs.default, resolve.defaults(list(x=pixdf),
                                              rest,
                                              list(labels=imnames, pch=".")))
      labels <- resolve.defaults(rest, list(labels=imnames))$labels
      colnames(pixdf) <- labels
    } else {
      do.call(hist.default,
              resolve.defaults(list(x=pixdf[,1L]),
                               rest,
                               list(main=paste("Histogram of", imnames[1L]),
                                    xlab=imnames[1L])))
    }
  }
  class(pixdf) <- unique(c("plotpairsim", class(pixdf)))
  attr(pixdf, "eps") <- eps
  return(invisible(pixdf))
}
