#
#   edges2triangles.R
#
#   $Revision: 1.9 $  $Date: 2013/05/21 09:40:26 $
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
                              nvert=max(iedge, jedge), check=TRUE) {
  if(check) {
    stopifnot(length(iedge) == length(jedge))
    stopifnot(length(iedge) == length(edgelength))
    stopifnot(all(iedge > 0))
    stopifnot(all(jedge > 0))
    if(!missing(nvert)) {
      stopifnot(all(iedge <= nvert))
      stopifnot(all(jedge <= nvert))
    }
  }
  # zero length data
  if(length(iedge) == 0)
    return(data.frame(i=integer(0),
                      j=integer(0),
                      k=integer(0),
                      diam=numeric(0)))

  # call C
  storage.mode(nvert) <- storage.mode(iedge) <- storage.mode(jedge) <- "integer"
  storage.mode(edgelength) <- "double"
  zz <- .Call("triDgraph",
              nv=nvert, iedge=iedge, jedge=jedge, edgelength=edgelength,
              PACKAGE="spatstat")
  df <- as.data.frame(zz)
  colnames(df) <- c("i", "j", "k", "diam")
  return(df)
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
  if(length(iedge) < 2) return(matrix(, nrow=0, ncol=3))
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

  
