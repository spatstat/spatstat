#
#  idw.R
#
#  Inverse-distance weighted smoothing
#
#  $Revision: 1.11 $ $Date: 2018/12/16 03:49:58 $

idw <- function(X, power=2, at=c("pixels", "points"), ..., se=FALSE) {
  stopifnot(is.ppp(X) && is.marked(X))
  at <- match.arg(at)
  marx <- marks(X)
  if(is.data.frame(marx)) {
    if((nc <- ncol(marx)) > 1) {
      ## multiple columns of marks - process one-by-one
      each <- vector(mode="list", length=nc)
      for(j in 1:nc) 
        each[[j]] <- idw(X %mark% marx[,j], power=power, at=at, ..., se=se)
      names(each) <- colnames(marx)
      ##
      if(!se) {
        ## estimates only
        switch(at,
               pixels = { out <- as.solist(each) },
               points = { out <- as.data.frame(each) }
               )
      } else {
        ## estimates and standard errors
        est <- lapply(each, getElement, name="estimate")
        SE  <- lapply(each, getElement, name="SE")
        switch(at,
               pixels = {
                 out <- list(estimate = as.solist(est),
                             SE       = as.solist(SE))
               },
               points = {
                 out <- list(estimate = as.data.frame(est),
                              SE      = as.data.frame(SE))
               })
      }
      return(out)
    } else 
      marx <- marx[,1L]
  }
  if(!is.numeric(marx))
    stop("Marks must be numeric")
  check.1.real(power)
  switch(at,
         pixels = {
           ## create grid
           W <- as.mask(as.owin(X), ...)
           dim <- W$dim
           npixels <- prod(dim)
           ## call C
           if(!se) {
             z <- .C("Cidw",
                     x = as.double(X$x),
                     y = as.double(X$y),
                     v = as.double(marx),
                     n = as.integer(npoints(X)),
                     xstart = as.double(W$xcol[1L]),
                     xstep  = as.double(W$xstep),
                     nx     = as.integer(dim[2L]),
                     ystart = as.double(W$yrow[1L]),
                     ystep  = as.double(W$ystep),
                     ny     = as.integer(dim[1L]),
                     power  = as.double(power),
                     num    = as.double(numeric(npixels)),
                     den    = as.double(numeric(npixels)),
                     rat    = as.double(numeric(npixels)),
                     PACKAGE = "spatstat")
             out <- as.im(matrix(z$rat, dim[1L], dim[2L]), W=W)
             out <- out[W, drop=FALSE]
           } else {
             z <- .C("Cidw2",
                     x = as.double(X$x),
                     y = as.double(X$y),
                     v = as.double(marx),
                     n = as.integer(npoints(X)),
                     xstart = as.double(W$xcol[1L]),
                     xstep  = as.double(W$xstep),
                     nx     = as.integer(dim[2L]),
                     ystart = as.double(W$yrow[1L]),
                     ystep  = as.double(W$ystep),
                     ny     = as.integer(dim[1L]),
                     power  = as.double(power),
                     num    = as.double(numeric(npixels)),
                     den    = as.double(numeric(npixels)),
                     rat    = as.double(numeric(npixels)),
                     mtwo   = as.double(numeric(npixels)),
                     wtwo   = as.double(numeric(npixels)),
                     PACKAGE = "spatstat")
             est <- as.im(matrix(z$rat, dim[1L], dim[2L]), W=W)
             est <- est[W, drop=FALSE]
             sumw  <- z$den
             sumw2 <- z$wtwo
             m2    <- z$mtwo
             varden <- sumw - sumw2/sumw
             varden[varden <= 0] <- NA
             SE <- sqrt(m2/varden)
             SE <- as.im(matrix(SE, dim[1L], dim[2L]), W=W)
             SE <- SE[W, drop=FALSE]
             out <- list(estimate=est, SE=SE)
           }
         },
         points={
           npts <- npoints(X)
           if(!se) {
             z <- .C("idwloo",
                     x = as.double(X$x),
                     y = as.double(X$y),
                     v = as.double(marx),
                     n = as.integer(npts),
                     power  = as.double(power),
                     num    = as.double(numeric(npts)),
                     den    = as.double(numeric(npts)),
                     rat    = as.double(numeric(npts)),
                     PACKAGE = "spatstat")
             out <- z$rat
           } else {
             z <- .C("idwloo2",
                     x = as.double(X$x),
                     y = as.double(X$y),
                     v = as.double(marx),
                     n = as.integer(npts),
                     power  = as.double(power),
                     num    = as.double(numeric(npts)),
                     den    = as.double(numeric(npts)),
                     rat    = as.double(numeric(npts)),
                     mtwo   = as.double(numeric(npts)),
                     wtwo   = as.double(numeric(npts)),
                     PACKAGE = "spatstat")
             est <- z$rat
             sumw  <- z$den
             sumw2 <- z$wtwo
             m2    <- z$mtwo
             varden <- sumw - sumw2/sumw
             varden[varden <= 0] <- NA
             SE <- sqrt(m2/varden)
             out <- list(estimate=est, SE=SE)
           }
         })
  return(out)
}
