#
# deldir.R
#
# Interface to deldir package
#
#  $Revision: 1.15 $ $Date: 2012/10/22 08:06:42 $
#

.spst.triEnv <- new.env()

assign("use.trigraf",  TRUE, envir=.spst.triEnv)
assign("use.trigrafS", TRUE, envir=.spst.triEnv)
assign("debug.delaunay", FALSE, envir=.spst.triEnv)

dirichlet <- function(X) {
  stopifnot(is.ppp(X))
  X <- unique(X, rule="deldir", warn=TRUE)
  w <- X$window
  dd <- deldir(X$x, X$y, rw=c(w$xrange,w$yrange))
  pp <- lapply(tile.list(dd), function(z) { owin(poly=z[c("x","y")]) })
  if(length(pp) == npoints(X))
    names(pp) <- seq_len(npoints(X))
  dir <- tess(tiles=pp, window=as.rectangle(w))
  if(w$type != "rectangle")
    dir <- intersect.tess(dir, w)
  return(dir)
}

delaunay <- function(X) {
  stopifnot(is.ppp(X))
  X <- unique(X, rule="deldir", warn=TRUE)
  nX <- npoints(X)
  if(nX < 3) return(NULL)
  w <- X$window
  dd <- deldir(X$x, X$y, rw=c(w$xrange, w$yrange))
  a <- dd$delsgs[,5]
  b <- dd$delsgs[,6]
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
            nt = as.integer(integer(1)),
            it = as.integer(integer(ne)),
            jt = as.integer(integer(ne)),
            kt = as.integer(integer(ne)),
            status = as.integer(integer(1)),
            PACKAGE="spatstat")
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
            nt = as.integer(integer(1)),
            it = as.integer(integer(ntmax)),
            jt = as.integer(integer(ntmax)),
            kt = as.integer(integer(ntmax)),
            status = as.integer(integer(1)),
            PACKAGE="spatstat")
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
  xtri <- matrix(x[tlist], nrow(tlist), 3)
  ytri <- matrix(y[tlist], nrow(tlist), 3)
  # ensure triangle vertices are in anticlockwise order
  ztri <- ytri - min(y)
  dx <- cbind(xtri[,2]-xtri[,1], xtri[,3]-xtri[,2], xtri[,1]-xtri[,3])
  zm <- cbind(ztri[,1]+ztri[,2], ztri[,2]+ztri[,3], ztri[,3]+ztri[,1])
  negareas <- apply(dx * zm, 1, sum)
  clockwise <- (negareas > 0)
  #
  if(any(clockwise)) {
    xc <- xtri[clockwise, , drop=FALSE]
    yc <- ytri[clockwise, , drop=FALSE]
    tc <- tlist[clockwise, , drop=FALSE]
    xtri[clockwise,]  <- xc[,c(1,3,2)]
    ytri[clockwise,]  <- yc[,c(1,3,2)]
    tlist[clockwise,] <- tc[, c(1,3,2)]
  }
  # At this point, triangle vertices are listed in anticlockwise order.
  # The same directed edge (i, j) cannot appear twice.
  # To weed out invalid triangles, check for such duplication
  triedges <- rbind(tlist[, c(1,2)],
                    tlist[, c(2,3)],
                    tlist[, c(3,1)])
  if(any(bad <- duplicated(triedges))) {
    badedges <- unique(triedges[bad, , drop=FALSE])
    ntri <- nrow(tlist)
    triid <- rep(seq_len(ntri), 3)
    illegal <- rep(FALSE, ntri)
    for(j in seq_len(nrow(badedges))) {
      from <- badedges[j, 1]
      to   <- badedges[j, 2]
      if(debug.delaunay)
        cat(paste("Suspect edge from vertex", from, "to vertex", to, "\n"))
      # find all triangles sharing this edge in this orientation
      sustri <- triid[(triedges[,1] == from) & (triedges[,2] == to)]
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

delaunay.distance <- function(X) {
  stopifnot(is.ppp(X))
  nX <- npoints(X)
  w <- as.owin(X)
  ok <- !duplicated(X, rule="deldir")
  Y <- X[ok] 
  nY <- npoints(Y)
  if(nY < 3) 
    return(matrix(Inf, nX, nX))
  dd <- deldir(Y$x, Y$y, rw=c(w$xrange,w$yrange))
  joins <- as.matrix(dd$delsgs[,5:6])
  joins <- rbind(joins, joins[,2:1])
  d <- matrix(-1L, nY, nY)
  diag(d) <- 0
  d[joins] <- 1
  adj <- matrix(FALSE, nY, nY)
  diag(adj) <- TRUE
  adj[joins] <- TRUE
  z <- .C("idist2dpath",
          nv = as.integer(nY),
          d = as.integer(d), 
          adj = as.integer(adj),
          dpath = as.integer(integer(nY * nY)),
          tol = as.integer(0),
          niter = as.integer(integer(1)), 
          status = as.integer(integer(1)),
          PACKAGE = "spatstat")
  if (z$status == -1)
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
