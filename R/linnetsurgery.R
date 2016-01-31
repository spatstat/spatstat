#'
#'   linnetsurgery.R
#'
#' Surgery on linear networks
#'
#' $Revision: 1.4 $  $Date: 2016/01/31 08:13:26 $
#'

insertVertices <- function(L, ...) {
  if(!inherits(L, c("lpp", "linnet")))
    stop("L should be a linear network (linnet) or point pattern (lpp)",
         call.=FALSE)
  if(haspoints <- is.lpp(L)) {
    X <- L
    L <- as.linnet(L)
    cooXnew <- cooXold <- coords(X) 
    segXold <- cooXold$seg 
    tpXold  <- cooXold$tp  
  }
  nsegL <- nsegments(L)
  ## validate new vertices
  V <- as.lpp(..., L=L)
  if(!identical(as.linnet(L, sparse=TRUE), as.linnet(V, sparse=TRUE)))
    stop("New vertices must lie on exactly the same network as L")
  if(npoints(V) == 0) {
    attr(L, "id") <- integer(0)
    if(!haspoints) {
      return(L)
    } else {
      X$domain <- L
      return(X)
    }
  }
  ## extract new vertex coordinates
  co <- coords(V)
  seg <- co$seg
  tp <- co$tp
  ## determine which segments will be split,
  ## and compute new serial numbers for the un-split segments
  splitsegments <- sort(unique(seg))
  notsplit <- rep(TRUE, nsegments(L))
  notsplit[splitsegments] <- FALSE
  segmap <- cumsum(notsplit)
  nunsplit <- sum(notsplit)
  ## existing vertices
  v <- L$vertices
  n <- npoints(v)
  ## initialise
  nadd <- 0
  vadd <- list(x=numeric(0), y=numeric(0))
  fromadd <- toadd <- id <- integer(0)
  ## split segments containing new vertices
  for(theseg in splitsegments) {
    ## find new vertices lying on segment 'theseg'
    i <- L$from[theseg]
    j <- L$to[theseg]
    those <- (seg == theseg)
    idthose <- which(those)
    ## order the new vertices along this segment
    tt <- tp[those]
    oo <- order(tt)
    tt <- tt[oo]
    idadd <- idthose[oo]
    ## make new vertices
    nnew <- length(tt)
    xnew <- with(v, x[i] + tt * diff(x[c(i,j)]))
    ynew <- with(v, y[i] + tt * diff(y[c(i,j)]))
    vnew <- list(x=xnew, y=ynew)
    ## make new edges
    kk <- n + nadd + (1:nnew)
    fromnew <- c(i, kk)
    tonew   <- c(kk, j)
    nnewseg <- nnew + 1
    ## add new vertices and edges to running total
    nadd <- nadd + nnew
    vadd <- concatxy(vadd, list(x=xnew, y=ynew))
    fromadd <- c(fromadd, fromnew)
    toadd <- c(toadd, tonew)
    id <- c(id, idadd)
    ## handle data points if any
    if(haspoints && any(relevant <- (segXold == theseg))) {
      tx <- tpXold[relevant]
      ttt <- c(0, tt, 1)
      m <- findInterval(tx, ttt, rightmost.closed=TRUE, all.inside=TRUE)
      t0 <- ttt[m]
      t1 <- ttt[m+1]
      tpXnew <- (tx - t0)/(t1-t0)
      tpXnew <- pmin(1, pmax(0, tpXnew))
      n0 <- nunsplit + length(fromadd) - nnewseg
      segXnew <- n0 + m
      cooXnew$seg[relevant] <- segXnew
      cooXnew$tp[relevant] <- tpXnew
    }
  }
  newfrom <- c(L$from[-splitsegments], fromadd)
  newto   <- c(L$to[-splitsegments], toadd)
  newv <- superimpose(v, vadd, check=FALSE)
  Lnew <- linnet(newv, edges=cbind(newfrom, newto),
                 sparse=identical(L$sparse, TRUE))
  newid <- integer(nadd)
  newid[id] <- n + 1:nadd
  attr(Lnew, "id") <- newid
  if(!haspoints)
    return(Lnew)
  ## adjust segment id for data points on segments that were not split
  Xnotsplit <- notsplit[segXold]
  cooXnew$seg[Xnotsplit] <- segmap[segXold[Xnotsplit]]
  Xnew <- lpp(cooXnew, Lnew)
  marks(Xnew) <- marks(X)
  attr(Xnew, "id") <- newid
  return(Xnew)
}
