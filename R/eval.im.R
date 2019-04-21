#
#     eval.im.R
#
#        eval.im()             Evaluate expressions involving images
#
#        compatible.im()       Check whether two images are compatible
#
#        harmonise.im()       Harmonise images
#        commonGrid()
#
#     $Revision: 1.54 $     $Date: 2019/04/21 11:42:27 $
#

eval.im <- local({

  eval.im <- function(expr, envir, harmonize=TRUE, warn=TRUE) {
    e <- as.expression(substitute(expr))
    ## get names of all variables in the expression
    varnames <- all.vars(e)
    allnames <- all.names(e, unique=TRUE)
    funnames <- allnames[!(allnames %in% varnames)]
    if(length(varnames) == 0)
      stop("No variables in this expression")
    ## get the values of the variables
    if(missing(envir)) {
      envir <- parent.frame() # WAS: sys.parent()
    } else if(is.list(envir)) {
      envir <- list2env(envir, parent=parent.frame())
    }
    vars <- mget(varnames, envir=envir, inherits=TRUE, ifnotfound=list(NULL))
    funs <- mget(funnames, envir=envir, inherits=TRUE, ifnotfound=list(NULL))
    ## WAS: vars <- lapply(as.list(varnames), get, envir=envir)
    ## WAS: funs <- lapply(as.list(funnames), get, envir=envir)
    ##
    ## find out which variables are images
    ims <- unlist(lapply(vars, is.im))
    if(!any(ims))
      stop("No images in this expression")
    images <- vars[ims]
    nimages <- length(images)
    ## test that the images are compatible
    if(!(do.call(compatible, unname(images)))) {
      whinge <- paste(if(nimages > 2) "some of" else NULL,
                      "the images",
                      commasep(sQuote(names(images))),
                      if(!harmonize) "are" else "were",
                      "not compatible")
      if(!harmonize) {
        stop(whinge, call.=FALSE)
      } else if(warn) {
        warning(whinge, call.=FALSE)
      }
      images <- do.call(harmonise.im, images)
    }
    ## trap a common error: using fv object as variable
    isfun <- unlist(lapply(vars, is.fv))
    if(any(isfun))
      stop("Cannot use objects of class fv as variables in eval.im")
    ## replace each image by its matrix of pixel values, and evaluate
    imagevalues <- lapply(images, getImValues)
    template <- images[[1L]]
    ## This bit has been repaired:
    vars[ims] <- imagevalues
    v <- eval(e, append(vars, funs))
    ##
    ## reshape, etc
    result <- im(v,
                 xcol=template$xcol, yrow=template$yrow,
                 xrange=template$xrange, yrange=template$yrange, 
                 unitname=unitname(template))
    return(result)
  }
  
  ## extract pixel values without destroying type information
  getImValues <- function(x) {
    v <- as.matrix(x)
    dim(v) <- NULL
    return(v)
  }

  eval.im
})

compatible.im <- function(A, B, ..., tol=1e-6) {
  verifyclass(A, "im")
  if(missing(B)) return(TRUE)
  verifyclass(B, "im")
  if(!all(A$dim == B$dim))
    return(FALSE)
  xdiscrep <- max(abs(A$xrange - B$xrange),
                 abs(A$xstep - B$xstep),
                 abs(A$xcol - B$xcol))
  ydiscrep <- max(abs(A$yrange - B$yrange),
                 abs(A$ystep - B$ystep),
                 abs(A$yrow - B$yrow))
  xok <- (xdiscrep < tol * min(A$xstep, B$xstep))
  yok <- (ydiscrep < tol * min(A$ystep, B$ystep))
  uok <- compatible.unitname(unitname(A), unitname(B))
  if(!(xok && yok && uok))
    return(FALSE)
  ## A and B are compatible
  if(length(list(...)) == 0)
    return(TRUE)
  ## recursion
  return(compatible.im(B, ..., tol=tol))
}

## force a list of images to be compatible

harmonize.im <- harmonise.im <- function(...) {
  argz <- list(...)
  n <- length(argz)
  if(n < 2) return(argz)
  result <- vector(mode="list", length=n)
  isim <- unlist(lapply(argz, is.im))
  if(!any(isim))
    stop("No images supplied")
  imgs <- argz[isim]
  ## if any windows are present, extract bounding box
  iswin <- unlist(lapply(argz, is.owin))
  bb0 <- if(!any(iswin)) NULL else do.call(boundingbox, unname(argz[iswin]))
  if(length(imgs) == 1L && is.null(bb0)) {
    ## only one 'true' image: use it as template.
    result[isim] <- imgs
    Wtemplate <- imgs[[1L]]
  } else {
    ## test for compatible units
    un <- lapply(imgs, unitname)
    uok <- unlist(lapply(un, compatible.unitname, y=un[[1L]]))
    if(!all(uok))
      stop("Images have incompatible units of length")
    ## find the image with the highest resolution
    xsteps <- unlist(lapply(imgs, getElement, name="xstep"))
    which.finest <- which.min(xsteps)
    finest <- imgs[[which.finest]]
    ## get the bounding box
    bb <- do.call(boundingbox, lapply(unname(imgs), as.rectangle))
    if(!is.null(bb0)) bb <- boundingbox(bb, bb0)
    ## determine new raster coordinates
    xcol <- prolongseq(finest$xcol, bb$xrange)
    yrow <- prolongseq(finest$yrow, bb$yrange)
    xy <- list(x=xcol, y=yrow)
    ## resample all images on new raster
    newimgs <- lapply(imgs, as.im, xy=xy)
    result[isim] <- newimgs
    Wtemplate <- newimgs[[which.finest]]
  }
  ## convert other data to images
  if(any(notim <- !isim)) 
    result[notim] <- lapply(argz[notim], as.im, W=as.mask(Wtemplate))
  names(result) <- names(argz)
  return(result)
}

## Return just the corresponding template window

commonGrid <- local({
  ## auxiliary function
  gettype <- function(x) {
    if(is.im(x) || is.mask(x)) "raster" else
    if(is.owin(x) || is.ppp(x) || is.psp(x)) "spatial" else
    "none"
  }

  commonGrid <- function(...) {
    argz <- list(...)
    type <- unlist(lapply(argz, gettype))
    israster <- (type == "raster")
    haswin   <- (type != "none")

    if(any(israster)) {
      ## Get raster data
      rasterlist <- argz[israster]
    } else {
      ## No existing raster data - apply default resolution
      if(!any(haswin))
        stop("No spatial data supplied")
      wins <- lapply(argz[haswin], as.owin)
      rasterlist <- lapply(wins, as.mask)
    }

    ## Find raster object with finest resolution
    if(length(rasterlist) == 1L) {
      ## only one raster object
      finest <- rasterlist[[1L]]
    } else {
      ## test for compatible units
      un <- lapply(rasterlist, unitname)
      uok <- unlist(lapply(un, compatible.unitname, y=un[[1L]]))
      if(!all(uok))
        stop("Objects have incompatible units of length")
      ## find the image/mask with the highest resolution
      xsteps <- unlist(lapply(rasterlist, getElement, name="xstep"))
      which.finest <- which.min(xsteps)
      finest <- rasterlist[[which.finest]]
    }
    ## determine the bounding box
    bb <- do.call(boundingbox, lapply(unname(argz[haswin]), as.rectangle))
    ## determine new raster coordinates
    xcol <- prolongseq(finest$xcol, bb$xrange)
    yrow <- prolongseq(finest$yrow, bb$yrange)
    xy <- list(x=xcol, y=yrow)
    ## generate template
    Wtemplate <- as.mask(bb, xy=xy)
    return(Wtemplate)
  }

  commonGrid
})

im.apply <- function(X, FUN, ..., fun.handles.na=FALSE, check=TRUE) {
  if(!inherits(X, "imlist")) {
    stopifnot(is.list(X))
    if(!all(sapply(X, is.im)))
      stop("All elements of X must be pixel images")
  }
  ## determine function to be applied
  fun <- if(is.character(FUN)) get(FUN, mode="function") else
         if(is.function(FUN)) FUN else stop("Unrecognised format for FUN")
  funcode <- match(list(fun),
                   list(base::sum,
                        base::mean,
                        base::mean.default,
                        stats::var,
                        stats::sd),
                   nomatch=0L)
  funtype <- c("general", "sum", "mean", "mean", "var", "sd")[funcode+1L]
  if(funcode != 0)
    na.rm <- resolve.1.default(list(na.rm=FALSE), list(...))
  ## ensure images are compatible
  if(check && !do.call(compatible, unname(X)))
    X <- do.call(harmonise.im, X)
  template <- X[[1L]]
  d <- dim(template)
  ## extract numerical values and convert to matrix with one column per image
  vals <- sapply(X, getElement, name="v")
  ## apply to all pixels ?
  full <- fun.handles.na || !anyNA(vals)
  if(!full) {
    ## NA present
    ok <- complete.cases(vals)
    if(!any(ok)) {
      ## empty result
      return(as.im(NA, W=template))
    }
    ## restrict to pixels where all data are non-NA
    vals <- vals[ok, , drop=FALSE]
  }
  n <- nrow(vals)
  ## calculate
  y <- switch(funtype,
              general = apply(vals, 1L, fun, ...),
              sum     = rowSums(vals,  na.rm=na.rm),
              mean    = rowMeans(vals, na.rm = na.rm),
              sd = ,
              var = {
                sumx  <- rowSums(vals,   na.rm = na.rm)
                sumx2 <- rowSums(vals^2, na.rm = na.rm)
                if(!anyNA(vals)) {
                  m <- ncol(vals)
                  v <- (sumx2 - sumx^2/m)/(m-1)
                } else {
                  m <- rowSums(!is.na(vals)) 
                  v <- ifelse(m < 2, NA, (sumx2 - sumx^2/m)/(m-1))
                }
                if(funtype == "var") v else sqrt(v)
              })
  if(funtype == "general" && length(y) != n)
    stop("FUN should yield one value per pixel")
  
  if(!full) {
    ## put the NA's back (preserving type of 'y')
    yfull <- rep(y[1L], prod(d))
    yfull[ok] <- y
    yfull[!ok] <- NA
    y <- yfull
  }
  
  ## pack up (preserving type of 'y')
  result <- im(y,
               xcol=template$xcol, yrow=template$yrow,
               xrange=template$xrange, yrange=template$yrange,
               unitname=template$unitname)
  return(result)
}
