#'
#'   linnetsurgery.R
#'
#' Surgery on linear networks and related objects
#'
#' $Revision: 1.31 $  $Date: 2020/04/01 04:44:15 $
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
  ## determine which segments will be split
  splitsegments <- sortunique(seg)
  notsplit <- rep(TRUE, nsegments(L))
  notsplit[splitsegments] <- FALSE
  ## un-split segments will be retained and listed first,
  ## followed by the new pieces.
  ## Compute new serial numbers for the un-split segments.
  segmap <- cumsum(notsplit) 
  nunsplit <- sum(notsplit)
  ## existing vertices
  v <- L$vertices
  n <- npoints(v)
  ## initialise lists of new vertices and edges 
  nadd <- 0
  vadd <- list(x=numeric(0), y=numeric(0))
  fromadd <- toadd <- id <- integer(0)
  ## initialise mapping from output segments to input segments
  comefrom <- which(notsplit) 
  ## split segments containing new vertices
  for(theseg in splitsegments) {
    ## find new vertices lying on segment 'theseg'
    those <- (seg == theseg)
    idthose <- which(those)
    ## order the new vertices along this segment
    tt <- tp[those]
    oo <- order(tt)
    tt <- tt[oo]
    idadd <- idthose[oo]
    ## make new vertices
    i <- L$from[theseg]
    j <- L$to[theseg]
    xnew <- with(v, x[i] + tt * diff(x[c(i,j)]))
    ynew <- with(v, y[i] + tt * diff(y[c(i,j)]))
    vnew <- list(x=xnew, y=ynew)
    nnew <- length(tt)
    ## replace edge i ~ j with edges i ~ k and k ~ j 
    kk <- n + nadd + seq_len(nnew)
    fromnew <- c(i, kk)
    tonew   <- c(kk, j)
    nnewseg <- nnew + 1L
    ## add new vertices and edges to running total
    nadd <- nadd + nnew
    vadd <- concatxy(vadd, vnew)
    fromadd <- c(fromadd, fromnew)
    toadd <- c(toadd, tonew)
    id <- c(id, idadd)
    comefrom <- c(comefrom, rep(theseg, nnewseg))
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
  edges <- cbind(newfrom, newto)
  reverse <- (newfrom > newto)
  if(anyrev <- any(reverse)) 
    edges[reverse, ] <- edges[reverse, c(2,1)]
  newv <- superimpose(v, vadd, check=FALSE)
  Lnew <- linnet(newv,
                 edges=edges,
                 sparse=identical(L$sparse, TRUE),
                 warn=FALSE)
  #' save information on provenance of line segments
  S <- Lnew$lines
  attr(S, "comefrom") <- comefrom
  if(length(comefrom) != nsegments(S))
    stop(paste("Internal error: length of provenance vector =",
               length(comefrom), "!=", nsegments(S), "= number of segments"),
         call.=FALSE)
  #' copy marks
  if(!is.null(marx <- marks(L$lines))) 
    marks(S) <- as.data.frame(marx)[comefrom, , drop=FALSE]
  Lnew$lines <- S
  #' save information identifying the new vertices in the new network
  newid <- integer(nadd)
  newid[id] <- n + 1:nadd
  attr(Lnew, "id") <- newid
  if(!haspoints)
    return(Lnew)
  ## adjust segment id for data points on segments that were not split
  Xnotsplit <- notsplit[segXold]
  cooXnew$seg[Xnotsplit] <- segmap[segXold[Xnotsplit]]
  ## adjust local coordinates if segment was reversed
  if(anyrev) {
    reversing <- reverse[cooXnew$seg]
    cooXnew$tp[reversing] <- 1 - cooXnew$tp[reversing]
  }
  Xnew <- lpp(cooXnew, Lnew)
  marks(Xnew) <- marks(X)
  attr(Xnew, "id") <- newid
  return(Xnew)
}

joinVertices <- function(L, from, to, marks=NULL) {
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
  Lnew <- linnet(vertices(L), edges=edges, sparse=L$sparse, warn=FALSE)
  #' assign marks only if provided
  if(!is.null(marks)) {
    marxnew <- (marks(L$lines) %mapp% marks) %msub% !duplicated(edges)
    if(!is.null(marxnew))
      marks(Lnew$lines) <- marxnew
  }
  if(!is.null(L$toler)) Lnew$toler <- L$toler
  if(!haspoints) return(Lnew)
  X <- lpp(Xdf, Lnew)
  return(X)
}

repairNetwork <- function(X) {
  if(!inherits(X, c("linnet", "lpp")))
    stop("X should be a linnet or lpp object", call.=FALSE)
  L <- as.linnet(X)
  V <- vertices(L)
  S    <- L$lines
  from <- L$from
  to   <- L$to
  ## check consistency between 'from,to' and 'lines'
  Sfrom <- endpoints.psp(S, "first")
  Sto   <- endpoints.psp(S, "second")
  fromV <- nncross(Sfrom, V, what="which")
  toV   <- nncross(Sto,   V, what="which")
  problem <- NULL
  if(length(from) != nsegments(S)) {
    problem <- "Network data implied different numbers of edges"
  } else if(any(to != toV) || any(from != fromV)) {
    #' find genuinely different vertices
    tol <- L$toler %orifnull% default.linnet.tolerance(S)
    clash <- (from != fromV)
    ss <- Sfrom[clash]
    vv <- V[from[clash]]
    d <- sqrt((ss$x - vv$x)^2 + (ss$y - vv$y)^2)
    bad <- (max(d) > tol)
    if(!bad) {
      clash <- (to != toV)
      ss <- Sto[clash]
      vv <- V[to[clash]]
      d <- sqrt((ss$x - vv$x)^2 + (ss$y - vv$y)^2)
      bad <- (max(d) > tol)
    }
    if(bad) 
      problem <- "Edge indices did not agree with segment endpoints"
  }
  if(!is.null(problem)) {
    if(is.marked(S)) {
      solution <- "edge indices were recomputed from line geometry"
      from <- L$from <- fromV
      to   <- L$to   <- toV
    } else {
      solution <- "lines were rebuilt from edge indices"
      xx <- V$x
      yy <- V$y
      L$lines <- psp(xx[from], yy[from], xx[to], yy[to], window=Window(V),
                     check=FALSE)
    }
    warning(paste0(problem, "; ", solution), call.=FALSE)
  } 

  reverse <- (from > to)
  if(any(reverse)) {
    newfrom <- ifelse(reverse, to, from)
    newto   <- ifelse(reverse, from, to)
    from <- L$from <- newfrom
    to   <- L$to   <- newto
    L$lines$ends[reverse,] <- L$lines$ends[reverse, c(3,4,1,2)]
    if(is.lpp(X)) { X$domain <- L } else { X <- L }
  }
  edgepairs <- cbind(from, to)
  retainedges <- !duplicated(as.data.frame(edgepairs)) & (from != to)
  keepall <- all(retainedges)
  if(is.lpp(X) && (!keepall || any(reverse))) {
    #' adjust segment coordinates
    cooX <- coords(X) # hyperframe, may include marks
    oldseg <- as.integer(unlist(cooX$seg))
    oldtp  <- as.numeric(unlist(cooX$tp))
    if(keepall) {
      newseg <- oldseg
    } else {
      segmap <- uniquemap(as.data.frame(edgepairs))
      newseg <- segmap[oldseg]
    }
    newtp  <- ifelse(reverse[oldseg], 1 - oldtp, oldtp)
    cooX$seg <- newseg
    cooX$tp <- newtp
    coords(X) <- cooX
  } 
  if(keepall)
    return(X)
  Y <- thinNetwork(X, retainedges=retainedges)
  return(Y)
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
  #' determine which edges/vertices are to be retained
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
  ## remove duplicate segments
  reverse <- (newfrom > newto)
  edgepairs <- cbind(ifelse(reverse, newto, newfrom),
                     ifelse(reverse, newfrom, newto))
  nontrivial <- (newfrom != newto) & !duplicated(edgepairs)
  edgepairs <- edgepairs[nontrivial,,drop=FALSE]
  reverse <- reverse[nontrivial]
  ## extract relevant subset of network
  Lsub <- linnet(Vsub, edges=edgepairs, sparse=sparse, warn=FALSE)
  ## tack on information about subset
  attr(Lsub, "retainvertices") <- retainvertices
  attr(Lsub, "retainedges") <- retainedges
  ## done?
  if(inherits(X, "linnet"))
    return(Lsub)
  ## X is an lpp object
  ## Find data points that lie on accepted segments
  dat <- X$data # hyperframe, may include marks
  ok <- retainedges[unlist(dat$seg)]
  dsub <- dat[ok, , drop=FALSE]
  ## compute new serial numbers for retained segments
  segmap <- cumsum(retainedges)
  oldseg <- as.integer(unlist(dsub$seg))
  dsub$seg <- newseg <- segmap[oldseg]
  ## adjust tp coordinate if segment endpoints were reversed
  if(any(revseg <- reverse[newseg])) {
    tp  <- as.numeric(unlist(dsub$tp))
    dsub$tp[revseg] <- 1 - tp[revseg]
  }
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

addVertices <- function(L, X, join=NULL, joinmarks=NULL) {
  if(!inherits(L, c("lpp", "linnet")))
    stop("L should be a linear network (linnet) or point pattern (lpp)",
         call.=FALSE)
  X <- as.ppp(X)
  if(haspoints <- is.lpp(L)) {
    Y <- L
    L <- as.linnet(L)
  }
  sparse <- L$sparse || is.null(L$dpath)
  V <- vertices(L)
  nV <- npoints(V)
  from <- L$from
  to   <- L$to
  ## new vertices
  nX <- npoints(X)
  Vplus <- superimpose(V, X, check=FALSE)
  nplus <- npoints(Vplus)
  iold <- seq_len(nV)
  inew <- nV + seq_len(nX)
  ## make new network
  Lplus <- L
  Lplus$vertices <- Vplus
  Lplus$window   <- Window(Vplus)
  Lplus$sparse   <- sparse
  ## 'lines', 'from', 'to', 'toler' are unchanged
  mplus <- sparseMatrix(i=c(from, to), j=c(to,from), x=TRUE,
                        dims=c(nplus, nplus))
  if(!sparse) mplus <- as.matrix(mplus)
  Lplus$m <- mplus
  if(!sparse) {
    dold <- L$dpath
    dnew <- matrix(Inf, nplus, nplus)
    diag(dnew) <- 0
    dnew[iold, iold] <- dold
    Lplus$dpath <- dnew
  }
  if(haspoints) {
    Y$domain <- Lplus # sufficient; point coordinates are still valid
  }
  out <- if(haspoints) Y else Lplus
  ## optionally join new vertices to existing network
  if(!is.null(join)) {
    if(is.numeric(join)) {
      check.nvector(join, nX, things="points of X")
      out <- joinVertices(out, inew, join, marks=joinmarks)
    } else if(is.character(join)) {
      join <- match.arg(join, c("vertices", "nearest"))
      switch(join,
             vertices={
               join <- nncross(X, V, what="which")
               out <- joinVertices(out, inew, join, marks=joinmarks)
             },
             nearest ={
               #' find nearest points on L
               p <- project2segment(X, as.psp(Lplus))
               locX <- data.frame(seg=p$mapXY, tp=p$tp)
               #' make them vertices
               out <- insertVertices(out, locX)
               #' join X to these new vertices
               joinid <- attr(out, "id")
               out <- joinVertices(out, inew, joinid, marks=joinmarks)
             })
    } else if(is.lpp(join)) {
      stopifnot(npoints(join) == npoints(X))
      out <- insertVertices(out, join)
      joinid <- attr(out, "id")
      out <- joinVertices(out, inew, joinid, marks=joinmarks)
    } else stop("Format of 'join' is not understood", call.=FALSE)
  }
  attr(out, "id") <- inew
  return(out)  
}
