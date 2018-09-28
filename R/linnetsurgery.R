#'
#'   linnetsurgery.R
#'
#' Surgery on linear networks and related objects
#'
#' $Revision: 1.13 $  $Date: 2018/06/12 02:55:59 $
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
  splitsegments <- sortunique(seg)
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
    vadd <- concatxy(vadd, vnew)
    fromadd <- c(fromadd, fromnew)
    toadd <- c(toadd, tonew)
    id <- c(id, idadd)
    ## handle data points if any
    if(haspoints && any(relevant <- (segXold == theseg))) {
      tx <- tpXold[relevant]
      ttt <- c(0, tt, 1)
      m <- findInterval(tx, ttt, rightmost.closed=TRUE, all.inside=TRUE)
      t0 <- ttt[m]
      t1 <- ttt[m+1L]
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

joinVertices <- function(L, from, to) {
  if(!inherits(L, c("lpp", "linnet")))
    stop("L should be a linear network (linnet) or point pattern (lpp)",
         call.=FALSE)
  if(haspoints <- is.lpp(L)) {
    X <- L
    L <- as.linnet(L)
    Xdf <- as.data.frame(X) 
  }
  if((missing(to) || is.null(to)) && !is.null(dim(from)) && ncol(from) == 2) {
    to   <- from[,2]
    from <- from[,1]
  } 
  newfrom <- as.integer(from)
  newto <- as.integer(to)
  edges <- cbind(c(L$from, newfrom), c(L$to, newto))
  Lnew <- linnet(vertices(L), edges=edges, sparse=L$sparse)
  if(!is.null(L$toler)) Lnew$toler <- L$toler
  if(!haspoints) return(Lnew)
  X <- lpp(Xdf, Lnew)
  return(X)
}

thinNetwork <- function(X, retainvertices, retainedges) {
  ## thin a network by retaining only the specified edges and/or vertices 
  if(!inherits(X, c("linnet", "lpp")))
    stop("X should be a linnet or lpp object", call.=FALSE)
  gotvert <- !missing(retainvertices)
  gotedge <- !missing(retainedges)
  if(!gotedge && !gotvert)
    return(X)
  L <- as.linnet(X)
  from <- L$from
  to   <- L$to
  V <- L$vertices
  sparse <- identical(L$sparse, TRUE)
  edgesFALSE <- logical(nsegments(L))
  verticesFALSE <- logical(npoints(V))
  if(!gotedge) {
    retainedges <- edgesFALSE
  } else if(!is.logical(retainedges)) {
    z <- edgesFALSE
    z[retainedges] <- TRUE
    retainedges <- z
  }
  if(!gotvert) {
    retainvertices <- verticesFALSE
  } else if(!is.logical(retainvertices)) {
    z <- verticesFALSE
    z[retainvertices] <- TRUE
    retainvertices <- z
  }
  if(gotvert && !gotedge) {
    ## retain all edges between retained vertices
    retainedges <- retainvertices[from] & retainvertices[to]
  } else if(gotedge) {
    ## retain vertices required for the retained edges
    retainvertices[from[retainedges]] <- TRUE
    retainvertices[to[retainedges]]   <- TRUE
  }
  ## assign new serial numbers to vertices, and recode
  Vsub <- V[retainvertices]
  newserial <- cumsum(retainvertices)
  newfrom <- newserial[from[retainedges]]
  newto   <- newserial[to[retainedges]]
  ## extract relevant subset of network
  Lsub <- linnet(Vsub, edges=cbind(newfrom, newto), sparse=sparse)
  ## tack on information about subset
  attr(Lsub, "retainvertices") <- retainvertices
  attr(Lsub, "retainedges") <- retainedges
  ## done?
  if(inherits(X, "linnet"))
    return(Lsub)
  ## X is an lpp object
  ## Find data points that lie on accepted segments
  dat <- X$data
  ok <- retainedges[unlist(dat$seg)]
  dsub <- dat[ok, , drop=FALSE]
  ## compute new serial numbers for retained segments
  segmap <- cumsum(retainedges)
  dsub$seg <- segmap[as.integer(unlist(dsub$seg))]
  # make new lpp object
  Y <- ppx(data=dsub, domain=Lsub, coord.type=as.character(X$ctype))
  class(Y) <- c("lpp", class(Y))
  ## tack on information about subset
  attr(Y, "retainpoints") <- ok
  return(Y)
}

validate.lpp.coords <- function(X, fatal=TRUE, context="") {
  ## check for mangled internal data
  proj <- project2segment(as.ppp(X), as.psp(as.linnet(X)))
  seg.claimed <- coords(X)$seg
  seg.mapped  <- proj$mapXY
  if(any(seg.claimed != seg.mapped)) {
    whinge <- paste("Incorrect segment id", context)
    if(fatal) stop(whinge, call.=FALSE) else warning(whinge, call.=FALSE)
    return(FALSE)
  }
  tp.claimed <- coords(X)$tp
  tp.mapped  <- proj$tp
  v <- max(abs(tp.claimed - tp.mapped))
  if(v > 0.01) {
    whinge <- paste("Incorrect 'tp' coordinate",
                    paren(paste("max discrepancy", v)),
                    context)
    if(fatal) stop(whinge, call.=FALSE) else warning(whinge, call.=FALSE)
    return(FALSE)
  }
  return(TRUE)
}
