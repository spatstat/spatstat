#'
#'   unstack.R
#'
#'   Methods for generic 'unstack'
#' 
#'   $Revision: 1.3 $  $Date: 2016/06/28 04:01:40 $

unstack.ppp <- unstack.psp <- unstack.lpp <- function(x, ...) {
  trap.extra.arguments(...)
  marx <- marks(x)
  d <- dim(marx)
  if(is.null(d)) return(solist(x))
  y <- rep(list(unmark(x)), d[2])
  for(j in seq_along(y))
    marks(y[[j]]) <- marx[,j,drop=FALSE]
  names(y) <- colnames(marx)
  return(as.solist(y))
}


unstack.msr <- function(x, ...) {
  trap.extra.arguments(...)
  d <- dim(x)
  if(is.null(d)) return(solist(x))
  smo <- attr(x, "smoothdensity")
  if(!inherits(smo, "imlist")) smo <- NULL
  nc <- d[2]
  y <- vector(mode="list", length=nc)
  for(j in seq_len(nc)) {
    xj <- x[,j,drop=FALSE]
    if(!is.null(smo)) attr(xj, "smoothdensity") <- smo[[j]]
    y[[j]] <- xj
  }
  names(y) <- colnames(x)
  return(as.solist(y))
}

unstackFilter <- function(x) {
  ## deal with a whole swag of classes that do not need to be unstacked
  nonvectorclasses <- c("im", "owin", "quad", "tess", 
                        "quadratcount", "quadrattest", 
                        "funxy", "distfun", "nnfun", 
                        "linnet", "linfun",
                        "influence.ppm", "leverage.ppm")
  y <- if(inherits(x, nonvectorclasses)) solist(x) else unstack(x)
  return(y)
}

unstack.solist <- function(x, ...) {
  trap.extra.arguments(...)
  as.solist(lapply(x, unstackFilter))
}

unstack.layered <- function(x, ...) {
  trap.extra.arguments(...)
  y <- lapply(x, unstackFilter)
  ny <- lengths(y)
  nx <- length(ny)
  if(all(ny == 1) || nx == 0) return(solist(x))
  pa <- layerplotargs(x)
  mm <- indexCartesian(ny)
  nz <- nrow(mm)
  z <- vector(mode="list", length=nz)
  nama <- lapply(y, names)
  for(i in seq_len(nz)) {
    ll <- mapply("[[", x=y, i=mm[i,], SIMPLIFY=FALSE)
    nam <- mapply("[", x=nama, i=mm[i,])
    nam <- nam[!sapply(nam, is.null)]
    names(z)[i] <- paste(nam, collapse=".")
    z[[i]] <- layered(LayerList=ll, plotargs=pa)
  }
  z <- as.solist(z)
  return(z)
}


  




  


  
