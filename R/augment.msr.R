#'
#'   augment.msr.R
#'
#'   Given a measure, compute a pixel image of the smoothed density
#'   and insert it in the object.
#'
#'   $Revision: 1.1 $  $Date: 2020/11/29 07:58:03 $


augment.msr <- function(x, ..., sigma, recompute=FALSE) {
  ## add a pixel image of the smoothed density component
  stopifnot(inherits(x, "msr"))
  if(!recompute && !is.null(attr(x, "smoothdensity"))) return(x)
  d <- ncol(as.matrix(x$val))
  xloc <- x$loc
  W <- as.owin(xloc)
  mt <- is.multitype(xloc)
  if(missing(sigma)) {
    sigma <- if(!mt) avenndist(xloc) else max(sapply(split(xloc), avenndist))
    if(sigma == 0) sigma <- max(bw.scott(xloc))/5
  }
  if(mt) {
    ## multitype case - split by type, extract smoothed part, then sum
    y <- lapply(split(x), augment.msr, sigma=sigma, ...)
    z <- lapply(y, attr, which="smoothdensity")
    if((nc <- ncol(x)) == 1) {
      ## scalar valued
      smo <- im.apply(z, sum)
      ## WAS:     z <- do.call(harmonise, unname(z))
      ##          smo <- Reduce("+", z)
    } else {
      ## vector valued
      smo <- vector(mode="list", length=nc)
      for(j in 1:nc) {
        zj <- lapply(z, "[[", i=j)
        smo[[j]] <- im.apply(zj, sum)
        ## WAS:    zj <- do.call(harmonise, unname(zj))
        ##         smo[[j]] <- Reduce("+", zj)
      }
      smo <- as.solist(smo)
    }
    attr(smo, "sigma") <- sigma
    attr(x, "smoothdensity") <- smo
    return(x)
  }
  ## Single-type 
  xdensity <- as.matrix(x$density)
  ## first weed out Inf, NA, NaN
  if(!all(ok <- complete.cases(xdensity))) 
    xdensity <- ok * xdensity
  ## smooth density unless constant
  ra <- apply(xdensity, 2, range)
  varble <- apply(as.matrix(ra), 2, diff) > sqrt(.Machine$double.eps)
  ##
  if(d == 1) {
    if(!varble) {
      smo <- as.im(mean(xdensity), W=W)
    } else {
      xmd <- xloc %mark% xdensity
      smo <- do.call(Smooth,
                     resolve.defaults(list(X=quote(xmd)),
                                      list(...),
                                      list(sigma=sigma)))
    }
  } else {
    smo <- vector(mode="list", length=d)
    names(smo) <- colnames(x)
    if(any(varble)) {
      xmdv <- xloc %mark% xdensity[,varble, drop=FALSE]
      smo[varble] <- do.call(Smooth,
                             resolve.defaults(list(X=quote(xmdv)),
                                              list(...),
                                              list(sigma=sigma)))
    }
    if(any(!varble)) 
      smo[!varble] <- solapply(apply(xdensity[, !varble, drop=FALSE], 2, mean),
                               as.im, W=W)
  }
  attr(smo, "sigma") <- sigma
  attr(x, "smoothdensity") <- smo
  return(x)
}

