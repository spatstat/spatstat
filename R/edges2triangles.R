#
#   edges2triangles.R
#
#   $Revision: 1.14 $  $Date: 2017/06/05 10:31:58 $
#

edges2triangles <- function(iedge, jedge, nvert=max(iedge, jedge),
                            ..., check=TRUE, friendly=rep(TRUE, nvert)) {
  usefriends <- !missing(friendly)
  if(check) {
    stopifnot(length(iedge) == length(jedge))
    stopifnot(all(iedge > 0))
    stopifnot(all(jedge > 0))
    if(!missing(nvert)) {
      stopifnot(all(iedge <= nvert))
      stopifnot(all(jedge <= nvert))
    }
    if(usefriends) {
      stopifnot(is.logical(friendly))
      stopifnot(length(friendly) == nvert)
      usefriends <- !all(friendly)
    }
  }
  # zero length data, or not enough to make triangles
  if(length(iedge) < 3) return(matrix(, nrow=0, ncol=3))
  # sort in increasing order of 'iedge'
  oi <- fave.order(iedge)
  iedge <- iedge[oi]
  jedge <- jedge[oi]
  # call C
  storage.mode(nvert) <- storage.mode(iedge) <- storage.mode(jedge) <- "integer"
  if(!usefriends) {
    zz <- .Call("triograph",
                nv=nvert, iedge=iedge, jedge=jedge,
                PACKAGE="spatstat")
  } else {
    fr <- as.logical(friendly)
    storage.mode(fr) <- "integer"
    zz <- .Call("trioxgraph",
                nv=nvert, iedge=iedge, jedge=jedge, friendly=fr,
                PACKAGE="spatstat")
  }
  mat <- as.matrix(as.data.frame(zz))
  return(mat)
}

# compute triangle diameters as well

trianglediameters <- function(iedge, jedge, edgelength, ..., 
                              nvert=max(iedge, jedge),
                              dmax=Inf, check=TRUE) {
  if(check) {
    stopifnot(length(iedge) == length(jedge))
    stopifnot(length(iedge) == length(edgelength))
    stopifnot(all(iedge > 0))
    stopifnot(all(jedge > 0))
    if(!missing(nvert)) {
      stopifnot(all(iedge <= nvert))
      stopifnot(all(jedge <= nvert))
    }
    if(is.finite(dmax)) check.1.real(dmax)
  }
  # zero length data
  if(length(iedge) == 0 || dmax < 0)
    return(data.frame(i=integer(0),
                      j=integer(0),
                      k=integer(0),
                      diam=numeric(0)))

  # call C
  storage.mode(nvert) <- storage.mode(iedge) <- storage.mode(jedge) <- "integer"
  storage.mode(edgelength) <- "double"
  if(is.infinite(dmax)) {
    zz <- .Call("triDgraph",
                nv=nvert, iedge=iedge, jedge=jedge, edgelength=edgelength,
                PACKAGE = "spatstat")
  } else {
    storage.mode(dmax) <- "double"
    zz <- .Call("triDRgraph",
                nv=nvert, iedge=iedge, jedge=jedge, edgelength=edgelength,
                dmax=dmax,
                PACKAGE = "spatstat")
  }    
  df <- as.data.frame(zz)
  colnames(df) <- c("i", "j", "k", "diam")
  return(df)
}

closetriples <- function(X, rmax) {
  a <- closepairs(X, rmax, what="ijd", twice=FALSE, neat=FALSE)
  tri <- trianglediameters(a$i, a$j, a$d, nvert=npoints(X), dmax=rmax)
  return(tri)
}

# extract 'vees', i.e. triples (i, j, k) where i ~ j and i ~ k

edges2vees <- function(iedge, jedge, nvert=max(iedge, jedge),
                            ..., check=TRUE) {
  if(check) {
    stopifnot(length(iedge) == length(jedge))
    stopifnot(all(iedge > 0))
    stopifnot(all(jedge > 0))
    if(!missing(nvert)) {
      stopifnot(all(iedge <= nvert))
      stopifnot(all(jedge <= nvert))
    }
  }
  # zero length data, or not enough to make vees
  if(length(iedge) < 2)
    return(data.frame(i=numeric(0),
                      j=numeric(0),
                      k=numeric(0)))
  # call 
  vees <- .Call("graphVees",
                nv = nvert,
                iedge = iedge,
                jedge = jedge,
                PACKAGE="spatstat")
  names(vees) <- c("i", "j", "k")
  vees <- as.data.frame(vees)
  return(vees)
}

  
