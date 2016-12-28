#
# deldir.R
#
# Interface to deldir package
#
#  $Revision: 1.27 $ $Date: 2015/10/21 09:06:57 $
#

.spst.triEnv <- new.env()

assign("use.trigraf",  TRUE, envir=.spst.triEnv)
assign("use.trigrafS", TRUE, envir=.spst.triEnv)
assign("debug.delaunay", FALSE, envir=.spst.triEnv)

dirichlet <- local({

  dirichlet <- function(X) {
    stopifnot(is.ppp(X))
    X <- unique(X, rule="deldir", warn=TRUE)
    w <- X$window
    dd <- safedeldir(X)
    if(is.null(dd)) return(NULL)
    pp <- lapply(tile.list(dd), df2poly)
    if(length(pp) == npoints(X))
      names(pp) <- seq_len(npoints(X))
    dir <- tess(tiles=pp, window=as.rectangle(w))
    if(w$type != "rectangle")
      dir <- intersect.tess(dir, w)
    return(dir)
  }

  df2poly <- function(z) { owin(poly=z[c("x","y")]) }
  
  dirichlet
})

delaunay <- function(X) {
  stopifnot(is.ppp(X))
  X <- unique(X, rule="deldir", warn=TRUE)
  nX <- npoints(X)
  if(nX < 3) return(NULL)
  w <- X$window
  dd <- safedeldir(X)
  if(is.null(dd)) return(NULL)
  a <- dd$delsgs[,5L]
  b <- dd$delsgs[,6L]
  use.trigraf  <- get("use.trigraf", envir=.spst.triEnv)
  use.trigrafS <- get("use.trigrafS", envir=.spst.triEnv)
  debug.delaunay <- get("debug.delaunay", envir=.spst.triEnv)
  if(use.trigrafS) {
    # first ensure a[] < b[]
    swap <- (a > b)
    if(any(swap)) {
      oldb <- b
      b[swap] <- a[swap]
      a[swap] <- oldb[swap]
    }
    # next ensure a is sorted
    o <- order(a, b)
    a <- a[o]
    b <- b[o]
    # 
    nv <- nX
    ne <- length(a)
    ntmax <- ne
    z <- .C("trigrafS",
            nv = as.integer(nv),
            ne = as.integer(ne),
            ie = as.integer(a),
            je = as.integer(b),
            ntmax = as.integer(ntmax),
            nt = as.integer(integer(1L)),
            it = as.integer(integer(ne)),
            jt = as.integer(integer(ne)),
            kt = as.integer(integer(ne)),
            status = as.integer(integer(1L)))
    if(z$status != 0)
      stop("Internal error: overflow in trigrafS")
    tlist <- with(z, cbind(it, jt, kt)[1:nt, ])
  } else if(use.trigraf) {
    nv <- nX
    ne <- length(a)
    ntmax <- ne
    z <- .C("trigraf",
            nv = as.integer(nv),
            ne = as.integer(ne),
            ie = as.integer(a),
            je = as.integer(b),
            ntmax = as.integer(ntmax),
            nt = as.integer(integer(1L)),
            it = as.integer(integer(ntmax)),
            jt = as.integer(integer(ntmax)),
            kt = as.integer(integer(ntmax)),
            status = as.integer(integer(1L)))
    if(z$status != 0)
      stop("Internal error: overflow in trigraf")
    tlist <- with(z, cbind(it, jt, kt)[1:nt, ])
  } else {
    tlist <- matrix(integer(0), 0, 3)
    for(i in seq_len(nX)) {
      # find all Delaunay neighbours of i 
      jj <- c(b[a==i], a[b==i])
      jj <- sort(unique(jj))
      # select those with a higher index than i
      jj <- jj[jj > i]
      # find pairs of neighbours which are Delaunay neighbours
      # (thus, triangles where the first numbered vertex is i)
      if(length(jj) > 0) 
        for(j in jj) {
          kk <- c(b[a == j], a[b == j])
          kk <- kk[(kk %in% jj) & (kk > j)]
          if(length(kk) > 0)
            for(k in kk) 
              # add (i,j,k) to list of triangles (i < j < k)
              tlist <- rbind(tlist, c(i, j, k))
        }
    }
  }
  # At this point, `tlist' contains all triangles formed by the Delaunay edges,
  # with vertices given in ascending order i < j < k in the 3 columns of tlist.
  # Some of these triangles may not belong to the Delaunay triangulation.
  # They will be weeded out later.
  
  # Assemble coordinates of triangles
  x <- X$x
  y <- X$y
  xtri <- matrix(x[tlist], nrow(tlist), 3L)
  ytri <- matrix(y[tlist], nrow(tlist), 3L)
  # ensure triangle vertices are in anticlockwise order
  ztri <- ytri - min(y)
  dx <- cbind(xtri[,2L]-xtri[,1L], xtri[,3L]-xtri[,2L], xtri[,1L]-xtri[,3L])
  zm <- cbind(ztri[,1L]+ztri[,2L], ztri[,2L]+ztri[,3L], ztri[,3L]+ztri[,1L])
  negareas <- apply(dx * zm, 1L, sum)
  clockwise <- (negareas > 0)
  #
  if(any(clockwise)) {
    xc <- xtri[clockwise, , drop=FALSE]
    yc <- ytri[clockwise, , drop=FALSE]
    tc <- tlist[clockwise, , drop=FALSE]
    xtri[clockwise,]  <- xc[,c(1L,3L,2L)]
    ytri[clockwise,]  <- yc[,c(1L,3L,2L)]
    tlist[clockwise,] <- tc[, c(1L,3L,2L)]
  }
  # At this point, triangle vertices are listed in anticlockwise order.
  # The same directed edge (i, j) cannot appear twice.
  # To weed out invalid triangles, check for such duplication
  triedges <- rbind(tlist[, c(1L,2L)],
                    tlist[, c(2L,3L)],
                    tlist[, c(3L,1L)])
  if(any(bad <- duplicated(triedges))) {
    badedges <- unique(triedges[bad, , drop=FALSE])
    ntri <- nrow(tlist)
    triid <- rep.int(seq_len(ntri), 3)
    illegal <- rep.int(FALSE, ntri)
    for(j in seq_len(nrow(badedges))) {
      from <- badedges[j, 1L]
      to   <- badedges[j, 2L]
      if(debug.delaunay)
        cat(paste("Suspect edge from vertex", from, "to vertex", to, "\n"))
      # find all triangles sharing this edge in this orientation
      sustri <- triid[(triedges[,1L] == from) & (triedges[,2L] == to)]
      if(debug.delaunay)
        cat(paste("\tInvestigating triangles", commasep(sustri), "\n"))
      # list all vertices associated with the suspect triangles
      susvert <- sort(unique(as.vector(tlist[sustri, ])))
      if(debug.delaunay)
        cat(paste("\tInvestigating vertices", commasep(susvert), "\n"))
      xsusvert <- x[susvert]
      ysusvert <- y[susvert]
      # take each triangle in turn and check whether it contains a data point
      for(k in sustri) {
        if(!illegal[k] &&
           any(inside.triangle(xsusvert, ysusvert, xtri[k,], ytri[k,]))) {
          if(debug.delaunay)
            cat(paste("Triangle", k, "is illegal\n"))
          illegal[k] <- TRUE
        }
      }
    }
    if(!any(illegal)) {
      if(debug.delaunay)
        cat("No illegal triangles found\n")
    } else {
      if(debug.delaunay)
        cat(paste("Removing", sum(illegal), "triangles\n"))
      tlist <- tlist[!illegal, , drop=FALSE]
      xtri  <- xtri[!illegal, , drop=FALSE]
      ytri  <- ytri[!illegal, , drop=FALSE]
    }
  }
  # make tile list
  tiles <- list()
  for(m in seq_len(nrow(tlist))) {
    p <- list(x=xtri[m,], y=ytri[m,])
    tiles[[m]] <- owin(poly=p, check=FALSE)
  }

  wc <- convexhull.xy(x, y)
  del <- tess(tiles=tiles, window=wc)
  if(w$type != "rectangle")
    del <- intersect.tess(del, w)
  return(del)
}

delaunayDistance <- function(X) {
  stopifnot(is.ppp(X))
  nX <- npoints(X)
  w <- as.owin(X)
  ok <- !duplicated(X, rule="deldir")
  Y <- X[ok] 
  nY <- npoints(Y)
  if(nY < 3) 
    return(matrix(Inf, nX, nX))
  dd <- deldir(Y$x, Y$y, rw=c(w$xrange,w$yrange))
  if(is.null(dd)) return(NULL)
  joins <- as.matrix(dd$delsgs[,5:6])
  joins <- rbind(joins, joins[,2:1])
  d <- matrix(-1L, nY, nY)
  diag(d) <- 0
  d[joins] <- 1
  adj <- matrix(FALSE, nY, nY)
  diag(adj) <- TRUE
  adj[joins] <- TRUE
  z <- .C("Idist2dpath",
          nv = as.integer(nY),
          d = as.integer(d), 
          adj = as.integer(adj),
          dpath = as.integer(integer(nY * nY)),
          tol = as.integer(0),
          niter = as.integer(integer(1L)), 
          status = as.integer(integer(1L)))
  if (z$status == -1L)
    warning(paste("graph connectivity algorithm did not converge after", 
                  z$niter, "iterations", "on", nY, "vertices and", 
                  sum(adj) - nY, "edges"))
  dpathY <- matrix(z$dpath, nY, nY)
  if(all(ok)) {
    dpathX <- dpathY
  } else {
    dpathX <- matrix(NA_integer_, nX, nX)
    dpathX[ok, ok] <- dpathY
  }
  return(dpathX)
}

safedeldir <- function(X) {
  rw <- with(X$window, c(xrange,yrange))
  dd <- try(deldir(X$x, X$y, rw=rw))
  if(!inherits(dd, "try-error") && inherits(dd, "deldir"))
    return(dd)
  warning("deldir failed; re-trying with slight perturbation of coordinates.",
          call.=FALSE)
  Y <- rjitter(X, mean(nndist(X))/100)
  dd <- try(deldir(Y$x, Y$y, rw=rw))
  if(!inherits(dd, "try-error") && inherits(dd, "deldir"))
    return(dd)
  warning("deldir failed even after perturbation of coordinates.", call.=FALSE)
  return(NULL)
}

dirichletVertices <- function(X) {
  DT <- tiles(dirichlet(X))
  xy <- do.call(concatxy, lapply(DT, vertices))
  Y <- unique(ppp(xy$x, xy$y, window=Window(X), check=FALSE))
  b <- bdist.points(Y)
  thresh <- diameter(Frame(X))/1000
  Y <- Y[b > thresh]
  return(Y)
}

dirichletAreas <- function(X) {
  stopifnot(is.ppp(X))
  X <- unmark(X)
  win <- Window(X)
  dup <- duplicated(X, rule="deldir")
  if((anydup <- any(dup))) {
    oldX <- X
    X <- X[!dup]
  }
  switch(win$type,
         rectangle = {
           rw <- c(win$xrange, win$yrange)
           dd <- deldir(X$x, X$y, dpl=NULL, rw=rw)
           w <- dd$summary[, 'dir.area']
         },
         polygonal = {
           w <- tile.areas(dirichlet(X))
         },
         mask = {
           #' Nearest data point to each pixel:
           tileid <- exactdt(X)$i
           #' Restrict to window (result is a vector - OK)
           tileid <- tileid[win$m]
           #' Count pixels in each tile
           id <- factor(tileid, levels=seq_len(X$n))
           counts <- table(id)
           #' Convert to digital area
           pixelarea <- win$xstep * win$ystep
           w <- pixelarea * as.numeric(counts)
         })
  if(!anydup)
    return(w)
  oldw <- numeric(npoints(oldX))
  oldw[!dup] <- w
  return(oldw)
}

delaunayNetwork <- function(X) {
  stopifnot(is.ppp(X))
  X <- unique(X, rule="deldir")
  nX <- npoints(X)
  if(nX == 0) return(NULL)
  if(nX == 1L) return(linnet(X, !diag(TRUE)))
  if(nX == 2L) return(linnet(X, !diag(c(TRUE,TRUE))))
  dd <- safedeldir(X)
  if(is.null(dd)) 
    return(NULL)
  joins <- as.matrix(dd$delsgs[, 5:6])
  return(linnet(X, edges=joins))
}

dirichletEdges <- function(X) {
  stopifnot(is.ppp(X))
  X <- unique(X, rule="deldir")
  nX <- npoints(X)
  W <- Window(X)
  if(nX < 2)
    return(edges(W))
  dd <- safedeldir(X)
  if(is.null(dd))
    return(edges(W))
  return(as.psp(dd$dirsgs[,1:4], window=W))
}

dirichletNetwork <- function(X, ...) as.linnet(dirichletEdges(X), ...)

## deprecated older names

delaunay.distance <- function(...) {
  .Deprecated("delaunayDistance", package="spatstat")
  delaunayDistance(...)
}

delaunay.network <- function(...) {
  .Deprecated("delaunayNetwork", package="spatstat")
  delaunayNetwork(...)
}

dirichlet.edges <- function(...) {
  .Deprecated("dirichletEdges", package="spatstat")
  dirichletEdges(...)
}

dirichlet.network <- function(...) {
  .Deprecated("dirichletNetwork", package="spatstat")
  dirichletNetwork(...)
}

dirichlet.vertices <- function(...) {
  .Deprecated("dirichletVertices", package="spatstat")
  dirichletVertices(...)
}
