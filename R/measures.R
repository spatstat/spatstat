#
#   measures.R
#
#  signed/vector valued measures with atomic and diffuse components
#
#  $Revision: 1.30 $  $Date: 2013/11/11 15:19:52 $
#
msr <- function(qscheme, discrete, density, check=TRUE) {
  if(!inherits(qscheme, "quad"))
    stop("qscheme should be a quadrature scheme")
  nquad <- n.quad(qscheme)
  U <- union.quad(qscheme)
  W <- w.quad(qscheme)
  Z <- is.data(qscheme)
  ndata <- sum(Z)
  # ensure conformable vectors/matrices
  if(is.vector(discrete) && is.vector(density)) {
    # handle constants
    if(length(discrete) == 1)
      discrete <- rep.int(discrete, ndata)
    if(length(density) == 1)
      density <- rep.int(density, nquad)
    # check lengths
    if(check) {
      check.nvector(discrete, ndata, things="data points", naok=TRUE)
      check.nvector(density,  nquad, things="quadrature points", naok=TRUE)
    }
    discretepad <- numeric(nquad)
    discretepad[Z] <- discrete
  } else {
    if(length(discrete) == 1 && is.matrix(density)) {
      # replicate constant 'discrete' component to matrix of correct size
      discrete <- matrix(discrete, ndata, ncol(density))
    } else if(length(density) == 1 && is.matrix(discrete)) {
      # replicate constant 'density' to matrix of correct size
      density <- matrix(density, nquad, ncol(discrete))
    } else {
      discrete <- as.matrix(discrete)
      density <- as.matrix(density)
    }
    if(check) {
      # check numbers of rows
      check.nmatrix(discrete, ndata, things="data points",
                    naok=TRUE, squarematrix=FALSE)
      check.nmatrix(density,  nquad, things="quadrature points",
                    naok=TRUE, squarematrix=FALSE)
    }
    nd <- ncol(discrete)
    nc <- ncol(density)
    if(nd != nc) {
      if(nd == 1) {
        # replicate columns of discrete component
        discrete <- matrix(rep.int(discrete, nc), ndata, nc)
        colnames(discrete) <- colnames(density)
      } else if(nc == 1) {
        # replicate columns of density component
        density <- matrix(rep.int(density, nd), nquad, nd)
        colnames(density) <- colnames(discrete)
      } else stop(paste("Incompatible numbers of columns in",
                        sQuote("discrete"), paren(nd), "and",
                        sQuote("density"), paren(nc)))
    }
    discretepad <- matrix(0, nquad, max(nd, nc))
    discretepad[Z, ] <- discrete
    colnames(discretepad) <- colnames(density)
  }

  #
  #
  # Discretised measure (value of measure for each quadrature tile)
  val <- discretepad + W * density
  if(is.matrix(density)) colnames(val) <- colnames(density)
  #
  out <- list(loc = U,
              val = val,
              atoms = Z,
              discrete = discretepad,
              density = density,
              wt = W)
  class(out) <- "msr"
  return(out)
}

# Translation table for usage of measures
#
#           e.g. res <- residuals(fit, ...)
#
#     OLD                               NEW           
#     res[ ]                       res$val[ ]       with(res, "increment")
#     attr(res, "atoms")           res$atoms        with(res, "is.atom")
#     attr(res, "discrete")        res$discrete     with(res, "discrete")
#     attr(res, "continuous")      res$density      with(res, "density")
#     w.quad(quad.ppm(fit))        res$wt           with(res, "qweights")
#     union.quad(quad.ppm(fit))    res$loc          with(res, "qlocations")
# .................................................

with.msr <- function(data, expr, ...) {
  stopifnot(inherits(data, "msr"))
  stopifnot(is.character(expr)) 
  y <- switch(expr,
              increment  = { data$val },
              is.atom    = { data$atoms },
              discrete   = { data$discrete },
              density    = { data$density },
              continuous = { data$density * data$wt },
              qweights   = { data$wt },
              qlocations = { data$loc },
              stop("Unrecognised option in entry.msr", call.=FALSE))
  return(y)
}

print.msr <- function(x, ...) {
  n <- npoints(x$loc)
  d <- ncol(as.matrix(x$val))
  descrip <- if(d == 1) "Scalar" else paste(d, "dimensional vector", sep="-")
  cat(paste(descrip, "-valued measure\n", sep=""))
  if(d > 1 && !is.null(cn <- colnames(x$val)))
    cat(paste("vector components:", commasep(sQuote(cn)), "\n"))
  cat(paste("Approximated by", n, "quadrature points\n"))
  print(as.owin(x$loc))
  cat(paste(sum(x$atoms), "atoms\n"))
  cat(paste("Total mass:\n"))
  if(d == 1) {
    cat(paste("discrete =", signif(sum(with(x, "discrete")), 5),
              "\tcontinuous =", signif(sum(with(x, "continuous")), 5),
              "\ttotal =", signif(sum(with(x, "increment")), 5), "\n"))
  } else {
    if(is.null(cn)) cn <- paste("component", 1:d)
    for(j in 1:d) {
      cat(paste(cn[j], ":\t",
                "discrete =", signif(sum(with(x, "discrete")[,j]), 5),
                "\tcontinuous =", signif(sum(with(x, "continuous")[,j]), 5),
                "\ttotal =", signif(sum(with(x, "increment")[,j]), 5), "\n"))
    }
  }
  return(invisible(NULL))
}

integral.msr <- function(x, ...) {
  stopifnot(inherits(x, "msr"))
  y <- with(x, "increment")
  if(is.matrix(y)) apply(y, 2, sum) else sum(y)
}
  
plot.msr <- function(x, ...) {
  xname <- short.deparse(substitute(x))
  d <- ncol(as.matrix(x$val))  
  if(d == 1) {
    # smooth the density unless it is flat
    if(diff(range(x$density)) > sqrt(.Machine$double.eps) ||
       "sigma" %in% names(list(...))) {
      sigma0 <- max(nndist(x$loc))
      smo <- do.call("Smooth",
                     resolve.defaults(list(X=x$loc %mark% x$density),
                                      list(...),
                                      list(sigma=sigma0)))
    } else {
      smo <- as.im(mean(x$density), W=as.owin(x$loc))
    }
    xtra <- unique(c(names(formals(plot.default)),
                     names(formals(image.default)),
                     "box"))
    do.call.matched("plot.im",
                    resolve.defaults(list(x=smo),
                                     list(...),
                                     list(main=xname)),
                    extrargs=xtra)
    xtra <- unique(c(names(formals(plot.owin)),
                     names(formals(points)),
                     names(formals(symbols))))
    xtra <- setdiff(xtra, "box")
    do.call.matched("plot.ppp",
                    resolve.defaults(list(x=x$loc %mark% x$discrete),
                                     list(add=TRUE),
                                     list(...)),
                    extrargs=xtra)
  } else {
    # split into a list of real-valued measures
    lis <- list()
    for(j in 1:d) 
      lis[[j]] <- x[,j]
    lis <- as.listof(lis)
    if(!is.null(cn <- colnames(x$val)))
      names(lis) <- cn
    do.call("plot.listof", resolve.defaults(list(lis),
                                            list(...),
                                            list(main=xname)))
  }
  return(invisible(NULL))  
}

"[.msr" <- function(x, i, j, ...) {
  valu  <- as.matrix(x$val)
  disc  <- as.matrix(x$discrete)
  dens  <- as.matrix(x$density)
  wt    <- x$wt
  atoms <- x$atoms
  #
  if(!missing(j)) {
    valu <- valu[, j]
    disc <- disc[, j]
    dens <- dens[, j]
  }
  loc <- x$loc
  if(!missing(i)) {
    # use [.ppp to identify which points are retained
    locn  <- loc %mark% seq_len(npoints(loc))
    loci  <- locn[i]
    loc   <- unmark(loci)
    id    <- marks(loci)
    # extract
    valu  <- valu[id, ]
    disc  <- disc[id, ]
    dens  <- dens[id, ]
    wt    <- wt[id]
    atoms <- atoms[id]
  }
  out <- list(loc=loc,
              val=valu,
              atoms=atoms,
              discrete=disc,
              density=dens,
              wt=wt)
  class(out) <- "msr"
  return(out)    
}

dim.msr <- function(x) { dim(as.matrix(x$val)) }

dimnames.msr <- function(x) { list(NULL, colnames(x$val)) }

smooth.msr <- function(X, ...) {
  message("smooth.msr will soon be deprecated: use the generic Smooth with a capital S")
#  .Deprecated("Smooth.msr", package="spatstat",
#     msg="smooth.msr is deprecated: use the generic Smooth with a capital S")
  Smooth(X, ...)
}

Smooth.msr <- function(X, ...) {
  verifyclass(X, "msr")
  loc <- X$loc
  val <- X$val
  d <- ncol(as.matrix(val))
  if(d == 1) {
    result <- density(loc, weights=val, ...)
  } else {
    result <- list()
    for(j in 1:d) 
      result[[j]] <- density(loc, weights=val[,j], ...)
    result <- as.listof(result)
    names(result) <- colnames(X)
  }
  return(result)
}
