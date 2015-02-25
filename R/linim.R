#
# linim.R
#
#  $Revision: 1.18 $   $Date: 2015/02/25 07:04:17 $
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
  cat("Image on linear network\n")
  print(attr(x, "L"))
  NextMethod("print")
}

plot.linim <- function(x, ..., style=c("colour", "width"), scale, adjust=1,
                       do.plot=TRUE) {
  xname <- short.deparse(substitute(x))
  style <- match.arg(style)
  # colour style: plot as pixel image
  if(style == "colour" || !do.plot)
    return(do.call("plot.im",
                   resolve.defaults(list(x),
                                    list(...),
                                    list(main=xname, do.plot=do.plot))))
  # width style
  L <- attr(x, "L")
  df <- attr(x, "df")
  Llines <- as.psp(L)
  # initialise plot
  W <- as.owin(L)
  bb <- do.call.matched("plot.owin",
                        resolve.defaults(list(x=W, type="n"),
                                         list(...), list(main=xname)),
                        extrargs="type")
  # resolve graphics parameters for polygons
  grafpar <- resolve.defaults(list(...), list(border=1, col=1))
  grafpar <- grafpar[names(grafpar) %in% names(formals(polygon))]
  # rescale values to a plottable range
  vr <- range(df$values)
  vr[1] <- min(0, vr[1])
  if(missing(scale)) {
    maxsize <- mean(distmap(Llines))/2
    scale <- maxsize/diff(vr)
  } 
  df$values <- adjust * scale * (df$values - vr[1])/2
  # split data by segment
  mapXY <- factor(df$mapXY, levels=seq_len(Llines$n))
  dfmap <- split(df, mapXY, drop=TRUE)
  # sort each segment's data by position along segment
  dfmap <- lapply(dfmap, function(z) { z[fave.order(z$tp), ] })
  # plot each segment's data
#  Lends <- Llines$ends
  Lperp <- angles.psp(Llines) + pi/2
  Lfrom <- L$from
  Lto   <- L$to
  Lvert <- L$vertices
  Ljoined  <- (vertexdegree(L) > 1)
  # precompute coordinates of dodecagon
  dodo <- disc(npoly=12)$bdry[[1]]
  #
  for(i in seq(length(dfmap))) {
    z <- dfmap[[i]]
    segid <- unique(z$mapXY)[1]
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
    vv <- c(vv[1],     vv, vv[length(vv)])
    rleft <- vv[1]
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
  return(invisible(bb))
}

as.im.linim <- function(X, ...) { as.im(X$Z, ...) }

as.linim <- function(X, ...) {
  UseMethod("as.linim")
}

as.linim.default <- function(X, L, ...) {
  stopifnot(inherits(L, "linnet"))
  Y <- as.im(X, W=as.rectangle(as.owin(L)), ...)
  Z <- as.im(as.mask.psp(as.psp(L), as.owin(Y)))
  Y <- eval.im(Z * Y)
  out <- linim(L, Y)
  return(out)
}

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
    envir <- sys.parent()
  } else if(is.list(envir)) {
    envir <- list2env(envir, parent=parent.frame())
  }
  vars <- lapply(as.list(varnames), function(x, e) get(x, envir=e), e=envir)
  names(vars) <- varnames
  funs <- lapply(as.list(funnames), function(x, e) get(x, envir=e), e=envir)
  names(funs) <- funnames
  # Find out which variables are (linear) images
  islinim <- unlist(lapply(vars, inherits, what="linim"))
  if(!any(islinim))
    stop("There are no linear images (class linim) in this expression")
  # ....................................
  # Evaluate the pixel values using eval.im
  # ....................................
  sc[[1]] <- as.name('eval.im')
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
      agree <- unlist(lapply(nets[-1], identical, y=nets[[1]]))
      if(!all(agree))
        stop(paste("Images do not refer to the same linear network"))
    }
    dfempty <- unlist(lapply(dframes, is.null))
    if(!any(dfempty)) {
      # replace each image variable by its data frame column of values
      vars[islinim] <- lapply(dframes, getElement, "values")
      # now evaluate expression
      Yvalues <- eval(e, append(vars, funs))
      # pack up
      dfY <- dframes[[1]]
      dfY$values <- Yvalues
    }
  }
  result <- linim(nets[[1]], Y, df=dfY)
  return(result)
}

as.linnet.linim <- function(X, ...) {
  attr(X, "L")
}

integral.linim <- function(f, domain=NULL, ...){
  verifyclass(f, "linim")
  L <- as.linnet(f)
  if(is.null(domain)) {
    mu <- mean(f)
    len <- volume(L)
  } else {
    mu <- mean(f[domain])
    len <- sum(lengths.psp(as.psp(L)[domain]))
  }
  return(mu * len)
}
