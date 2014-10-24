# 
# linnet.R
#    
#    Linear networks
#
#    $Revision: 1.25 $    $Date: 2014/10/24 00:22:30 $
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

linnet <- function(vertices, m, edges) {
  if(missing(m) && missing(edges))
    stop("specify either m or edges")
  if(!missing(m) && !missing(edges))
    stop("do not specify both m and edges")
  # validate inputs
  stopifnot(is.ppp(vertices))
  if((nv <- npoints(vertices)) <= 1) {
    m <- matrix(FALSE, nv, nv)
  } else if(!missing(m)) {
    # check logical matrix
    stopifnot(is.matrix(m) && is.logical(m) && isSymmetric(m))
    if(nrow(m) != vertices$n)
      stop("dimensions of matrix m do not match number of vertices")
  } else {
    # check (from, to) pairs
    stopifnot(is.matrix(edges) && ncol(edges) == 2)
    if(any((edges %% 1) != 0))
      stop("Entries of edges list should be integers")
    np <- npoints(vertices)
    if(any(edges > np))
      stop("index out-of-bounds in edges list")
    # convert to adjacency matrix
    m <- matrix(FALSE, np, np)
    m[edges] <- TRUE
    m <- m | t(m)
  }
  # create line segments
  rowm <- row(m)
  colm <- col(m)
  uptri <- m & (rowm < colm)
  from <- as.vector(rowm[uptri])
  to   <- as.vector(colm[uptri])
  xx   <- vertices$x
  yy   <- vertices$y
  lines <- psp(xx[from], yy[from], xx[to], yy[to], window=vertices$window,
               check=FALSE)
  # compute matrix of distances between adjacent vertices
  n <- nrow(m)
  d <- matrix(Inf, n, n)
  diag(d) <- 0
  d[m] <- pairdist(vertices)[m]
  # now compute shortest-path distances between each pair of vertices
  dpath <- dist2dpath(d)
  if(any(is.infinite(dpath)))
    warning("Network is not connected")
  # pack up
  out <- list(vertices=vertices, m=m, lines=lines, from=from, to=to,
              dpath=dpath, window=vertices$window)
  class(out) <- c("linnet", class(out))
  # pre-compute circumradius
  out$circumradius <- circumradius(out)
  return(out)  
}

print.linnet <- function(x, ...) {
  nv <- x$vertices$n
  nl <- x$lines$n
  ne <- sum(x$m)/2
  splat("Linear network with",
        nv, ngettext(nv, "vertex,", "vertices,"), 
        nl, ngettext(nl, "line", "lines"),
        "and",
        ne, ngettext(ne, "edge", "edges"))
  return(invisible(NULL))
}

summary.linnet <- function(object, ...) {
  print(object, ...)
  unitinfo <- summary(unitname(object))
  cat(paste("Total length",
            sum(lengths.psp(object$lines)),
            unitinfo$plural, unitinfo$explain, "\n"))
  print(as.owin(object))
  return(invisible(NULL))
}

plot.linnet <- function(x, ..., main=NULL, add=FALSE,
                        vertices=FALSE, window=FALSE) {
  if(is.null(main))
    main <- short.deparse(substitute(x))
  stopifnot(inherits(x, "linnet"))
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
  plot(lines, ..., add=TRUE, main=main)
  if(vertices)
    plot(x$vertices, add=TRUE)
  return(invisible(NULL))
}

as.psp.linnet <- function(x, ..., fatal=TRUE) {
  verifyclass(x, "linnet", fatal=fatal)
  return(x$lines)
}

as.owin.linnet <- function(W, ...) {
  return(as.owin(as.psp(W)))
}

as.linnet <- function(X, ...) {
  UseMethod("as.linnet")
}

as.linnet.linnet <- function(X, ...) { X }

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
  max(0, x$dpath)
}

volume.linnet <- function(x) {
  sum(lengths.psp(x$lines))
}

circumradius.linnet <- function(x, ...) {
  stopifnot(inherits(x, "linnet"))
  cr <- x$circumradius
  if(!is.null(cr))
    return(cr)
  dpath <- x$dpath
  if(nrow(dpath) <= 1)
    return(max(0,dpath))
  from  <- x$from
  to    <- x$to
  lines <- x$lines
  nseg  <- lines$n
  leng  <- lengths.psp(lines)
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
  Y$dpath        <- f * X$dpath
  Y$circumradius <- f * X$circumradius
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
    Y$dpath        <- scal * X$dpath
    Y$circumradius <- scal * X$circumradius
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

"[.linnet" <- function(x, i, ...) {
  if(!is.owin(i))
    stop("In [.linnet: the index i should be a window", call.=FALSE)
  # Find vertices that lie inside 'i'
  okvert <- inside.owin(x$vertices, w=i)
  # find segments whose endpoints both lie in 'i'
  okedge <- okvert[x$from] & okvert[x$to]
  # assign new serial numbers to vertices, and recode 
  newserial <- cumsum(okvert)
  newfrom <- newserial[x$from[okedge]]
  newto   <- newserial[x$to[okedge]]
  # make new linear network
  xnew <- linnet(x$vertices[i], edges=cbind(newfrom, newto))
  return(xnew)
}
