# 
# linnet.R
#    
#    Linear networks
#
#    $Revision: 1.19 $    $Date: 2011/03/21 09:08:24 $
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
  if(!missing(m)) {
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
  cat(paste("Linear network with",
            x$vertices$n, "vertices,",
            x$lines$n, "lines and",
            sum(x$m/2), "edges\n"))
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
  max(x$dpath)
}

circumradius <- function(x) {
  stopifnot(inherits(x, "linnet"))
  cr <- x$circumradius
  if(!is.null(cr))
    return(cr)
  dpath <- x$dpath
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
