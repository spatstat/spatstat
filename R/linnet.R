# 
# linnet.R
#    
#    Linear networks
#
#    $Revision: 1.61 $    $Date: 2017/02/07 08:12:05 $
#
# An object of class 'linnet' defines a linear network.
# It includes the following components
#
#        vertices     (ppp)      vertices of network
#
#        m            (matrix)   adjacency matrix
#
#        lines        (psp)      edges of network
#
#        dpath        (matrix)   matrix of shortest path distances
#                                between each pair of vertices
#
#        from, to     (vectors)  map from edges to vertices.
#                                The endpoints of the i-th segment lines[i]
#                                are vertices[from[i]] and vertices[to[i]]
#
#
#  FUNCTIONS PROVIDED:
#       linnet        creates an object of class "linnet" from data
#       print.linnet  print an object of class "linnet"
#       plot.linnet   plot an object of class "linnet"
#

# Make an object of class "linnet" from the minimal data

linnet <- function(vertices, m, edges, sparse=FALSE, warn=TRUE) {
  if(missing(m) && missing(edges))
    stop("specify either m or edges")
  if(!missing(m) && !missing(edges))
    stop("do not specify both m and edges")
  # validate inputs
  stopifnot(is.ppp(vertices))
  nv <- npoints(vertices)
  if(nv <= 1) {
    m <- matrix(FALSE, nv, nv)
  } else if(!missing(m)) {
    # check logical matrix or logical sparse matrix
    if(!is.matrix(m) && !inherits(m, c("lgCMatrix", "lgTMatrix")))
      stop("m should be a matrix or sparse matrix")
    stopifnot(is.logical(m) && isSymmetric(m))
    if(nrow(m) != vertices$n)
      stop("dimensions of matrix m do not match number of vertices")
    if(any(diag(m))) {
      warning("diagonal entries of the matrix m should not be TRUE; ignored")
      diag(m) <- FALSE
    }
    sparse <- !is.matrix(m)
    ## determine 'from' and 'to' vectors
    ij <- which(m, arr.ind=TRUE)
    ij <- ij[ ij[,1L] < ij[,2L], , drop=FALSE]
    from <- ij[,1L]
    to   <- ij[,2L]
  } else {
    # check (from, to) pairs
    stopifnot(is.matrix(edges) && ncol(edges) == 2)
    if(any((edges %% 1) != 0))
      stop("Entries of edges list should be integers")
    if(any(self <- (edges[,1L] == edges[,2L]))) {
      warning("edge list should not join a vertex to itself; ignored")
      edges <- edges[!self, , drop=FALSE]
    }
    np <- npoints(vertices)
    if(any(edges > np))
      stop("index out-of-bounds in edges list")
    from <- edges[,1L]
    to   <- edges[,2L]
    # convert to adjacency matrix
    if(!sparse) {
      m <- matrix(FALSE, np, np)
      m[edges] <- TRUE
    } else 
      m <- sparseMatrix(i=from, j=to, x=TRUE, dims=c(np, np))
    m <- m | t(m)
  }
  # create line segments
  xx   <- vertices$x
  yy   <- vertices$y
  lines <- psp(xx[from], yy[from], xx[to], yy[to], window=vertices$window,
               check=FALSE)
  # tolerance
  toler <- default.linnet.tolerance(lines)
  ## pack up
  out <- list(vertices=vertices, m=m, lines=lines, from=from, to=to,
              sparse=sparse, window=vertices$window,
              toler=toler)
  class(out) <- c("linnet", class(out))
  ## finish ?
  if(sparse)
    return(out)
  # compute matrix of distances between adjacent vertices
  n <- nrow(m)
  d <- matrix(Inf, n, n)
  diag(d) <- 0
  d[m] <- pairdist(vertices)[m]
  ## now compute shortest-path distances between each pair of vertices
  out$dpath <- dpath <- dist2dpath(d)
  if(warn && any(is.infinite(dpath)))
    warning("Network is not connected", call.=FALSE)
  # pre-compute bounding radius 
  out$boundingradius <- boundingradius(out)
  return(out)  
}

print.linnet <- function(x, ...) {
  nv <- x$vertices$n
  nl <- x$lines$n
  splat("Linear network with",
        nv, ngettext(nv, "vertex", "vertices"), 
        "and",
        nl, ngettext(nl, "line", "lines"))
  print(as.owin(x), prefix="Enclosing window: ")
  return(invisible(NULL))
}

summary.linnet <- function(object, ...) {
  deg <- vertexdegree(object)
  sparse <- object$sparse %orifnull% is.null(object$dpath)
  result <- list(nvert = object$vertices$n,
                 nline = object$lines$n,
                 nedge = sum(deg)/2,
                 unitinfo = summary(unitname(object)),
                 totlength = sum(lengths.psp(object$lines)),
                 maxdegree = max(deg),
                 win = as.owin(object),
                 sparse = sparse)
  if(!sparse) {
    result$diam <- diameter(object)
    result$boundrad <- boundingradius(object)
  }
  result$toler <- object$toler
  class(result) <- c("summary.linnet", class(result))
  result
}

print.summary.linnet <- function(x, ...) {
  dig <- getOption('digits')
  with(x, {
    splat("Linear network with",
          nvert, ngettext(nvert, "vertex", "vertices"), 
          "and",
          nline, ngettext(nline, "line", "lines"))
    splat("Total length", signif(totlength, dig), 
          unitinfo$plural, unitinfo$explain)
    splat("Maximum vertex degree:", maxdegree)
    if(sparse) splat("[Sparse matrix representation]") else {
      splat("Diameter:", signif(diam, dig), unitinfo$plural)
      splat("Bounding radius:", signif(boundrad, dig), unitinfo$plural)
    }
    if(!is.null(x$toler))
      splat("Numerical tolerance:", signif(x$toler, dig), unitinfo$plural)
    print(win, prefix="Enclosing window: ")
  })
  return(invisible(NULL))
}

plot.linnet <- function(x, ..., main=NULL, add=FALSE,
                        vertices=FALSE, window=FALSE,
                        do.plot=TRUE) {
  if(is.null(main))
    main <- short.deparse(substitute(x))
  stopifnot(inherits(x, "linnet"))
  bb <- Frame(x)
  if(!do.plot) return(invisible(bb))
  lines <- as.psp(x)
  if(!add) {
    # initialise new plot
    w <- as.owin(lines)
    if(window)
      plot(w, ..., main=main)
    else
      plot(w, ..., main=main, type="n")
  }
  # plot segments and (optionally) vertices
  do.call(plot,
          resolve.defaults(list(x=lines,
                                show.all=FALSE, add=TRUE,
                                main=main),
                           list(...)))
  if(vertices)
    plot(x$vertices, add=TRUE)
  return(invisible(bb))
}

as.psp.linnet <- function(x, ..., fatal=TRUE) {
  verifyclass(x, "linnet", fatal=fatal)
  return(x$lines)
}

vertices.linnet <- function(w) {
  verifyclass(w, "linnet")
  return(w$vertices)
}

nvertices.linnet <- function(x, ...) {
  verifyclass(x, "linnet")
  return(x$vertices$n)
}

nsegments.linnet <- function(x) {
  return(x$lines$n)
}

Window.linnet <- function(X, ...) {
  return(X$window)
}

"Window<-.linnet" <- function(X, ..., check=TRUE, value) {
  if(check) {
    X <- X[value]
  } else {
    X$window <- value
    X$lines$window <- value
    X$vertices$window <- value
  }
  return(X)
}

as.owin.linnet <- function(W, ...) {
  return(Window(W))
}

as.linnet <- function(X, ...) {
  UseMethod("as.linnet")
}

as.linnet.linnet <- function(X, ..., sparse) {
  if(missing(sparse)) return(X)
  if(is.null(X$sparse)) X$sparse <- is.null(X$dpath)
  if(sparse && !(X$sparse)) {
    # delete distance matrix
    X$dpath <- NULL
    # convert adjacency matrix to sparse matrix
    X$m <- as(X$m, "sparseMatrix")
    X$sparse <- TRUE
  } else if(!sparse && X$sparse) {
    # convert adjacency to matrix
    X$m <- m <- as.matrix(X$m)
    edges <- which(m, arr.ind=TRUE)
    from <- edges[,1L]
    to   <- edges[,2L]
    # compute distances to one-step neighbours
    n <- nrow(m)
    d <- matrix(Inf, n, n)
    diag(d) <- 0
    coo <- coords(vertices(X))
    d[edges] <- sqrt(rowSums((coo[from, 1:2] - coo[to, 1:2])^2))
    # compute shortest path distance matrix
    X$dpath <- dist2dpath(d)
    # compute bounding radius
    X$boundingradius <- boundingradius(X)
    X$sparse <- FALSE
  } else if(!sparse) {
    # possibly update internals
    X$boundingradius <- boundingradius(X)
  }
  # possibly update internals
  X$circumradius <- NULL
  X$toler <- default.linnet.tolerance(X)
  return(X)
}

as.linnet.psp <- local({
  
  as.linnet.psp <- function(X, ..., eps, sparse=FALSE) {
    X <- selfcut.psp(X)
    V <- unique(endpoints.psp(X))
    if(missing(eps) || is.null(eps)) {
      eps <- sqrt(.Machine$double.eps) * diameter(Frame(X))
    } else {
      check.1.real(eps)
      stopifnot(eps >= 0)
    }
    if(eps > 0 && minnndist(V) <= eps) {
      gV <- marks(connected(V, eps))
      xy <- split(coords(V), gV)
      mxy <- lapply(xy, centro)
      V <- do.call(superimpose, append(unname(mxy), list(W=Window(X))))
    }
    first  <- endpoints.psp(X, "first")
    second <- endpoints.psp(X, "second")
    from <- nncross(first, V, what="which")
    to   <- nncross(second, V, what="which")
    nontrivial <- (from != to)
    join <- cbind(from, to)[nontrivial, , drop=FALSE]
    result <- linnet(V, edges=join, sparse=sparse)
    if(is.marked(X)) marks(result$lines) <- marks(X[nontrivial])
    return(result)
  }

  centro <- function(X) as.list(apply(X, 2, mean))
  
  as.linnet.psp
})


unitname.linnet <- function(x) {
  unitname(x$window)
}

"unitname<-.linnet" <- function(x, value) {
  w <- x$window
  v <- x$vertices
  l <- x$lines
  unitname(w) <- unitname(v) <- unitname(l) <- value
  x$window <- w
  x$vertices <- v
  x$lines <- l
  return(x)
}

diameter.linnet <- function(x) {
  stopifnot(inherits(x, "linnet"))
  dpath <- x$dpath
  if(is.null(dpath)) return(NULL) else return(max(0, dpath))
}

volume.linnet <- function(x) {
  sum(lengths.psp(x$lines))
}

vertexdegree <- function(x) {
  verifyclass(x, "linnet")
  return(rowSums(x$m))
}

circumradius.linnet <- function(x, ...) {
  .Deprecated("boundingradius.linnet")
  boundingradius.linnet(x, ...)
}

boundingradius.linnet <- function(x, ...) {
  stopifnot(inherits(x, "linnet"))
  cr <- x$boundingradius %orifnull% x$circumradius
  if(!is.null(cr))
    return(cr)
  dpath <- x$dpath
  if(is.null(dpath)) return(NULL)
  if(any(is.infinite(dpath))) return(Inf)
  if(nrow(dpath) <= 1)
    return(max(0,dpath))
  from  <- x$from
  to    <- x$to
  lines <- x$lines
  nseg  <- lines$n
  leng  <- lengths.psp(lines)
  if(spatstat.options("Clinearradius")) {
    fromC <- from - 1L
    toC   <- to - 1L
    nv <- npoints(vertices(x))
    huge <- sum(leng)
    z <- .C("linearradius",
            ns = as.integer(nseg),
            from = as.integer(fromC),
            to = as.integer(toC),
            lengths = as.double(leng),
            nv = as.integer(nv), 
            dpath = as.double(dpath), 
            huge = as.double(huge), 
            result = as.double(numeric(1)),
            PACKAGE = "spatstat")
    return(z$result)
  }
  sA <- sB <- matrix(Inf, nseg, nseg)
  for(i in 1:nseg) {
    # endpoints of segment i
    A <- from[i]
    B <- to[i]
    AB <- leng[i]
    sA[i,i] <- sB[i,i] <- AB/2
    for(j in (1:nseg)[-i]) {
    # endpoints of segment j
      C <- from[j]
      D <- to[j]
      CD <- leng[j]
      AC <- dpath[A,C]
      AD <- dpath[A,D]
      BC <- dpath[B,C]
      BD <- dpath[B,D]
      # max dist from A to any point in segment j
      sA[i,j] <- if(AD > AC + CD) AC + CD else
                if(AC > AD + CD) AD + CD else
                (AC + AD + CD)/2
      # max dist from B to any point in segment j
      sB[i,j] <- if(BD > BC + CD) BC + CD else
                if(BC > BD + CD) BD + CD else
                (BC + BD + CD)/2
    }
  }
  # max dist from each A to any point in another segment
  mA <- apply(sA, 1, max)
  # max dist from each B to any point in another segment
  mB <- apply(sB, 1, max)
  # min of these
  min(mA, mB)
}



####################################################
# affine transformations
####################################################

scalardilate.linnet <- function(X, f, ...) {
  trap.extra.arguments(..., .Context="In scalardilate(X,f)")
  check.1.real(f, "In scalardilate(X,f)")
  stopifnot(is.finite(f) && f > 0)
  Y <- X
  Y$vertices     <- scalardilate(X$vertices, f=f)
  Y$lines        <- scalardilate(X$lines, f=f)
  Y$window       <- scalardilate(X$window, f=f)
  if(!is.null(X$dpath)) {
    Y$dpath        <- f * X$dpath
    Y$boundingradius <- f * (X$boundingradius %orifnull% X$circumradius)
    Y$circumradius <- NULL
  }
  if(!is.null(X$toler))
    X$toler <- makeLinnetTolerance(f * X$toler)
  return(Y)
}

affine.linnet <- function(X,  mat=diag(c(1,1)), vec=c(0,0), ...) {
  verifyclass(X, "linnet")
  if(length(unique(eigen(mat)$values)) == 1) {
    # transformation is an isometry
    scal <- sqrt(abs(det(mat)))
    Y <- X
    Y$vertices     <- affine(X$vertices, mat=mat, vec=vec, ...)
    Y$lines        <- affine(X$lines,    mat=mat, vec=vec, ...)
    Y$window       <- affine(X$window,   mat=mat, vec=vec, ...)
    if(!is.null(X$dpath)) {
      Y$dpath        <- scal * X$dpath
      Y$boundingradius <- scal * (X$boundingradius %orifnull% X$circumradius)
      X$circumradius <- NULL
    }
    if(!is.null(Y$toler))
      Y$toler <- makeLinnetTolerance(scal * Y$toler)
  } else {
    # general case
    vertices <- affine(X$vertices, mat=mat, vec=vec, ...)
    Y <- linnet(vertices, edges=cbind(X$from, X$to))
  }
  return(Y)
}

shift.linnet <- function(X, vec=c(0,0), ..., origin=NULL) {
  verifyclass(X, "linnet")
  Y <- X
  Y$window  <- W <- shift(X$window, vec=vec, ..., origin=origin)
  v <- getlastshift(W)
  Y$vertices <- shift(X$vertices, vec=v, ...)
  Y$lines    <- shift(X$lines, vec=v, ...)
  # tack on shift vector
  attr(Y, "lastshift") <- v
  return(Y)
}

rotate.linnet <- function(X, angle=pi/2, ..., centre=NULL) {
  verifyclass(X, "linnet")
  if(!is.null(centre)) {
    X <- shift(X, origin=centre)
    negorigin <- getlastshift(X)
  } else negorigin <- NULL
  Y <- X
  Y$vertices <- rotate(X$vertices, angle=angle, ...)
  Y$lines    <- rotate(X$lines, angle=angle, ...)
  Y$window   <- rotate(X$window, angle=angle, ...)
  if(!is.null(negorigin))
    Y <- shift(Y, -negorigin)
  return(Y)
}

rescale.linnet <- function(X, s, unitname) {
  if(missing(unitname)) unitname <- NULL
  if(missing(s) || is.null(s)) s <- 1/unitname(X)$multiplier
  Y <- scalardilate(X, f=1/s)
  unitname(Y) <- rescale(unitname(X), s, unitname)
  return(Y)
}

"[.linnet" <- function(x, i, ..., snip=TRUE) {
  if(!is.owin(i))
    stop("In [.linnet: the index i should be a window", call.=FALSE)
  w <- i
  ## Find vertices that lie inside window
  vertinside <- inside.owin(x$vertices, w=w)
  from <- x$from
  to   <- x$to
  if(snip) {
    ## For efficiency, first restrict network to relevant segments.
    ## Find segments EITHER OF whose endpoints lie in 'w'
    okedge <- vertinside[from] | vertinside[to]
    ## extract relevant subset of network graph
    x <- thinNetwork(x, retainedges=okedge)
    ## Now add vertices at crossing points with boundary of 'w'
    b <- crossing.psp(as.psp(x), edges(w))
    x <- insertVertices(x, unique(b))
    boundarypoints <- attr(x, "id")
    ## update data
    from <- x$from
    to   <- x$to
    vertinside <- inside.owin(x$vertices, w=w)
    vertinside[boundarypoints] <- TRUE
  }
  ## find segments whose endpoints BOTH lie in 'w'
  edgeinside <- vertinside[from] & vertinside[to]
  ## extract relevant subset of network
  xnew <- thinNetwork(x, retainedges=edgeinside)
  ## adjust window efficiently
  Window(xnew, check=FALSE) <- w
  return(xnew)
}

#
# interactive plot for linnet objects
#

iplot.linnet <- function(x, ..., xname) {
  if(missing(xname))
    xname <- short.deparse(substitute(x))
  if(!inherits(x, "linnet"))
    stop("x should be a linnet object")
  ## predigest
  v <- vertices(x)
  deg <- vertexdegree(x)
  dv <- textstring(v, txt=paste(deg))
  y <- layered(lines=as.psp(x),
               vertices=v,
               degree=dv)
  iplot(y, ..., xname=xname, visible=c(TRUE, FALSE, FALSE))
}

pixellate.linnet <- function(x, ...) {
  pixellate(as.psp(x), ...)
}

connected.linnet <- function(X, ..., what=c("labels", "components")) {
  verifyclass(X, "linnet")
  what <- match.arg(what)
  nv <- npoints(vertices(X))
  ie <- X$from - 1L
  je   <- X$to   - 1L
  ne <- length(ie)
  zz <- .C("cocoGraph",
           nv = as.integer(nv),
           ne = as.integer(ne), 
           ie = as.integer(ie),
           je = as.integer(je),
           label = as.integer(integer(nv)), 
           status = as.integer(integer(1L)),
           PACKAGE = "spatstat")
  if (zz$status != 0) 
    stop("Internal error: connected.linnet did not converge")
  lab <- zz$label + 1L
  lab <- as.integer(factor(lab))
  lab <- factor(lab)
  if(what == "labels")
    return(lab)
  nets <- list()
  subsets <- split(seq_len(nv), lab)
  for(i in seq_along(subsets)) 
    nets[[i]] <- thinNetwork(X, retainvertices=subsets[[i]])
  return(nets)
}

crossing.linnet <- function(X, Y) {
  X <- as.linnet(X)
  if(!inherits(Y, c("linnet", "infline", "psp")))
    stop("L should be an object of class psp, linnet or infline", call.=FALSE)
  ## convert infinite lines to segments
  if(inherits(Y, "linnet")) Y <- as.psp(Y)
  if(inherits(Y, "infline")) {
    Y <- clip.infline(Y, Frame(X))
    id <- marks(Y)
    lev <- levels(id)
  } else {
    id <- lev <- seq_len(nsegments(Y))
  }
  ## extract segments of network
  S <- as.psp(X)
  ## find crossing points
  SY <- crossing.psp(S, Y, fatal=FALSE, details=TRUE)
  if(is.null(SY) || npoints(SY) == 0)
    return(lpp(L=X))
  SY <- as.data.frame(SY)
  Z <- with(as.data.frame(SY),
            as.lpp(x=x, y=y, seg=iA, tp=tA, L=X,
                   marks=factor(id[as.integer(jB)], levels=lev)))
  return(Z)
}
