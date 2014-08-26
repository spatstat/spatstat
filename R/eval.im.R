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
#     $Revision: 1.32 $     $Date: 2014/07/29 09:04:19 $
#

eval.im <- function(expr, envir, harmonize=TRUE) {
  e <- as.expression(substitute(expr))
  # get names of all variables in the expression
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
  # find out which variables are images
  ims <- unlist(lapply(vars, is.im))
  if(!any(ims))
    stop("No images in this expression")
  images <- vars[ims]
  nimages <- length(images)
  # test that the images are compatible
  if(!(ok <- do.call("compatible", unname(images)))) {
    whinge <- paste(if(nimages > 2) "some of" else NULL,
                    "the images",
                    commasep(sQuote(names(images))),
                    if(!harmonize) "are" else "were",
                    "not compatible")
    if(!harmonize) {
      stop(whinge, call.=FALSE)
    } else {
      warning(whinge, call.=FALSE)
      images <- do.call("harmonise.im", images)
    }
  }
  # replace each image by its matrix of pixel values, and evaluate
  getvalues <- function(x) {
    v <- as.matrix(x)
    dim(v) <- NULL
    return(v)
  }
  imagevalues <- lapply(images, getvalues)
  template <- images[[1]]
  # This bit has been repaired:
  vars[ims] <- imagevalues
  v <- eval(e, append(vars, funs))
  #
  # reshape, etc
  result <- im(v,
               xcol=template$xcol, yrow=template$yrow,
               xrange=template$xrange, yrange=template$yrange, 
               unitname=unitname(template))
  return(result)
}
  
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
  uok <- compatible.units(unitname(A), unitname(B))
  if(!(xok && yok && uok))
    return(FALSE)
  # A and B are compatible
  if(length(list(...)) == 0)
    return(TRUE)
  # recursion
  return(compatible.im(B, ..., tol=tol))
}

# force a list of images to be compatible

harmonize.im <- harmonise.im <- function(...) {
  argz <- list(...)
  n <- length(argz)
  if(n < 2) return(argz)
  result <- vector(mode="list", length=n)
  isim <- unlist(lapply(argz, is.im))
  if(!any(isim))
    stop("No images supplied")
  imgs <- argz[isim]
  # if any windows are present, extract bounding box
  iswin <- unlist(lapply(argz, is.owin))
  bb0 <- if(!any(iswin)) NULL else do.call("boundingbox", unname(argz[iswin]))
  if(length(imgs) == 1 && is.null(bb0)) {
    # only one 'true' image: use it as template.
    result[isim] <- imgs
    Wtemplate <- imgs[[1]]
  } else {
    # test for compatible units
    un <- lapply(imgs, unitname)
    uok <- unlist(lapply(un, compatible.units, y=un[[1]]))
    if(!all(uok))
      stop("Images have incompatible units of length")
    # find the image with the highest resolution
    xsteps <- unlist(lapply(imgs, function(a) { a$xstep }))
    which.finest <- which.min(xsteps)
    finest <- imgs[[which.finest]]
    # get the bounding box
    bb <- do.call("boundingbox", lapply(unname(imgs), as.rectangle))
    if(!is.null(bb0)) bb <- boundingbox(bb, bb0)
    # determine new raster coordinates
    xcol <- prolongseq(finest$xcol, bb$xrange)
    yrow <- prolongseq(finest$yrow, bb$yrange)
    xy <- list(x=xcol, y=yrow)
    # resample all images on new raster
    newimgs <- lapply(imgs, as.im, xy=xy)
    result[isim] <- newimgs
    Wtemplate <- newimgs[[which.finest]]
  }
  # convert other data to images
  if(any(notim <- !isim)) 
    result[notim] <- lapply(argz[notim], as.im, W=as.mask(Wtemplate))
  names(result) <- names(argz)
  return(result)
}

# Return just the corresponding template window

commonGrid <- local({
  # auxiliary function
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
      # Get raster data
      rasterlist <- argz[israster]
    } else {
      # No existing raster data - apply default resolution
      if(!any(haswin))
        stop("No spatial data supplied")
      wins <- lapply(argz[haswin], as.owin)
      rasterlist <- lapply(wins, as.mask)
    }

    # Find raster object with finest resolution
    if(length(rasterlist) == 1) {
      # only one raster object
      finest <- rasterlist[[1]]
    } else {
      # test for compatible units
      un <- lapply(rasterlist, unitname)
      uok <- unlist(lapply(un, compatible.units, y=un[[1]]))
      if(!all(uok))
        stop("Objects have incompatible units of length")
      # find the image/mask with the highest resolution
      xsteps <- unlist(lapply(rasterlist, function(a) { a$xstep }))
      which.finest <- which.min(xsteps)
      finest <- rasterlist[[which.finest]]
    }
    # determine the bounding box
    bb <- do.call("boundingbox", lapply(unname(argz[haswin]), as.rectangle))
    # determine new raster coordinates
    xcol <- prolongseq(finest$xcol, bb$xrange)
    yrow <- prolongseq(finest$yrow, bb$yrange)
    xy <- list(x=xcol, y=yrow)
    # generate template
    Wtemplate <- as.mask(bb, xy=xy)
    return(Wtemplate)
  }

  commonGrid
})

im.apply <- function(X, FUN, ...) {
  stopifnot(is.list(X))
  if(!all(unlist(lapply(X, is.im))))
    stop("All elements of imlist must be pixel images")
  fun <- if(is.character(FUN)) get(FUN) else
         if(is.function(FUN)) FUN else stop("Unrecognised format for FUN")
  ## ensure images are compatible
  X <- do.call(harmonise.im, X)
  template <- X[[1]]
  ## extract numerical values and convert to matrix with one column per image
  vlist <- lapply(X, function(z) as.vector(as.matrix(z)))
  vals <- matrix(unlist(vlist), ncol=length(X))
  colnames(vals) <- names(X)
  ok <- complete.cases(vals)
  if(!any(ok)) {
    ## empty result
    return(as.im(NA, W=template))
  }
  ## apply function
  resultok <- apply(vals[ok,, drop=FALSE], 1, fun, ...)
  if(length(resultok) != sum(ok))
    stop("FUN should yield one value per pixel")
  ## pack up, with attention to type of data
  d <- dim(template)
  resultmat <- matrix(resultok[1], d[1], d[2])
  resultmat[ok] <- resultok
  resultmat[!ok] <- NA
  result <- as.im(resultmat, W=X[[1]])
  if(is.factor(resultok)) levels(result) <- levels(resultok)
  return(result)
}
