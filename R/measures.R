#
#   measures.R
#
#  signed/vector valued measures with atomic and diffuse components
#
#  $Revision: 1.47 $  $Date: 2014/07/28 08:26:13 $
#
msr <- function(qscheme, discrete, density, check=TRUE) {
  if(!inherits(qscheme, "quad"))
    stop("qscheme should be a quadrature scheme")
  nquad <- n.quad(qscheme)
  U <- union.quad(qscheme)
  wt <- w.quad(qscheme)
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
  val <- discretepad + wt * density
  if(is.matrix(density)) colnames(val) <- colnames(density)
  #
  out <- list(loc = U,
              val = val,
              atoms = Z,
              discrete = discretepad,
              density = density,
              wt = wt)
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
  splat(paste0(if(d == 1) "Scalar" else paste0(d, "-dimensional vector"),
               "-valued measure"))
  if(d > 1 && !is.null(cn <- colnames(x$val)))
    splat("vector components:", commasep(sQuote(cn)))
  if(waxlyrical("extras")) {
    splat("Approximated by", n, "quadrature points")
    print(as.owin(x$loc))
    splat(sum(x$atoms), "atoms")
    splat("Total mass:")
    if(d == 1) {
      splat("discrete =", signif(sum(with(x, "discrete")), 5),
            "  continuous =", signif(sum(with(x, "continuous")), 5),
            "  total =", signif(sum(with(x, "increment")), 5))
    } else {
      if(is.null(cn)) cn <- paste("component", 1:d)
      for(j in 1:d) {
        splat(paste0(cn[j], ":\t"),
              "discrete =", signif(sum(with(x, "discrete")[,j]), 5),
              "  continuous =", signif(sum(with(x, "continuous")[,j]), 5),
              "  total =", signif(sum(with(x, "increment")[,j]), 5))
      }
    }
  }
  return(invisible(NULL))
}

integral.msr <- function(x, ...) {
  stopifnot(inherits(x, "msr"))
  y <- with(x, "increment")
  if(is.matrix(y)) apply(y, 2, sum) else sum(y)
}

augment.msr <- function(x, ..., sigma) {
  ## add a pixel image of the smoothed density component
  stopifnot(inherits(x, "msr"))
  d <- ncol(as.matrix(x$val))
  xloc <- x$loc
  W <- as.owin(xloc)
  if(missing(sigma)) sigma <- maxnndist(xloc)
  ## smooth density unless constant
  xdensity <- as.matrix(x$density)
  ra <- apply(xdensity, 2, range)
  varble <- apply(as.matrix(ra), 2, diff) > sqrt(.Machine$double.eps)
  ##
  if(d == 1) {
    smo <- if(!varble) as.im(mean(xdensity), W=W) else
           do.call("Smooth",
                   resolve.defaults(list(X=xloc %mark% xdensity),
                                    list(...),
                                    list(sigma=sigma)))
  } else {
    smo <- vector(mode="list", length=d)
    names(smo) <- colnames(x)
    if(any(varble)) 
      smo[varble] <-
        do.call("Smooth",
                resolve.defaults(list(X=xloc %mark% xdensity[,varble]),
                                 list(...),
                                 list(sigma=sigma)))
    if(any(!varble)) 
      smo[!varble] <- lapply(apply(xdensity[, !varble], 2, mean),
                             function(z, W) as.im(z, W=W),
                             W=W)
  }
  attr(x, "smoothdensity") <- smo
  return(x)
}

plot.msr <- function(x, ..., add=FALSE,
                     how=c("image", "contour", "imagecontour"),
                     do.plot=TRUE) {
  xname <- short.deparse(substitute(x))
  how <- match.arg(how)
  d <- ncol(as.matrix(x$val))
  xloc <- x$loc
  if(is.null(smo <- attr(x, "smoothdensity"))) {
    x <- augment.msr(x, ...)
    smo <- attr(x, "smoothdensity")
  }
  if(d > 1) {
    ## split into a list of real-valued measures
    lis <- list()
    for(j in 1:d) {
      xj <- x[,j]
      attr(xj, "smoothdensity") <- smo[[j]]
      lis[[j]] <- xj
    }
    lis <- as.listof(lis)
    if(!is.null(cn <- colnames(x$val)))
      names(lis) <- cn
    result <- do.call("plot.listof", resolve.defaults(list(lis),
                                                      list(...),
                                                      list(how=how,
                                                           main=xname)))
    return(invisible(result))
  }
  ## 
  xatomic <- (x$loc %mark% x$discrete)[x$atoms]
  xtra.im <- unique(c(names(formals(plot.default)),
                      names(formals(image.default)),
                      "box"))
  xtra.pp <- unique(c(names(formals(plot.owin)),
                      "pch", "cex", "lty", "lwd",
                      names(formals(symbols)),
                      "etch"))
  xtra.pp <- setdiff(xtra.pp, "box")
  ##
  do.image <-  how %in% c("image", "imagecontour")
  do.contour <-  how %in% c("contour", "imagecontour")
  ## allocate space for plot and legend using do.plot=FALSE mechanism
  pdata <- do.call.matched("plot.ppp",
                           resolve.defaults(list(x=xatomic,
                                                 do.plot=FALSE),
                                            list(...)),
                           extrargs=xtra.pp)
  result <- pdata
  bb <- attr(pdata, "bbox")
  if(do.image) {
    idata <- do.call.matched("plot.im",
                             resolve.defaults(list(x=smo, do.plot=FALSE),
                                              list(...)),
                             extrargs=xtra.im)
    result <- idata
    bb <- boundingbox(bb, attr(idata, "bbox"))
  }
  ##
  attr(result, "bbox") <- bb
  ##
  if(do.plot) {
    if(!add) 
      ## initialise plot
      do.call.matched(plot.owin,
                      resolve.defaults(list(x=bb, type="n", main="   "),
                                       list(...)))
    ## display density
    if(do.image) 
      do.call.matched("plot.im",
                      resolve.defaults(list(x=smo, add=TRUE),
                                       list(...),
                                       list(main=xname, show.all=TRUE)),
                      extrargs=xtra.im)
    if(do.contour) 
      do.call.matched("contour.im",
                      resolve.defaults(list(x=smo, add=TRUE),
                                       list(...),
                                       list(main=xname,
                                            axes=FALSE, show.all=!do.image)),
                      extrargs=c("zlim", "labels", "labcex",
                        ## DO NOT ALLOW 'col' 
                        "drawlabels", "method", "vfont", "lty", "lwd"))
    ## display atoms
    do.call.matched("plot.ppp",
                    resolve.defaults(list(x=xatomic, add=TRUE, main=""),
                                     list(...)),
                    extrargs=xtra.pp)
  }
  return(invisible(result))
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
  .Deprecated("Smooth.msr", package="spatstat",
     msg="smooth.msr is deprecated: use the generic Smooth with a capital S")
  Smooth(X, ...)
}

Smooth.msr <- function(X, ..., drop=TRUE) {
  verifyclass(X, "msr")
  loc <- X$loc
  val <- X$val
  result <- density(loc, weights=val, ...)
  if(!drop && is.im(result))
    result <- listof(result)
  return(result)
}

as.owin.msr <- function(W, ..., fatal=TRUE) {
  as.owin(W$loc, ..., fatal=fatal)
}

domain.msr <- Window.msr <- function(X, ...) { as.owin(X) } 

shift.msr <- function(X,  ...) {
  X$loc <- Xloc <- shift(X$loc, ...)
  putlastshift(X, getlastshift(Xloc))
}

as.layered.msr <- function(X) {
  nc <- ncol(X)
  if(nc == 0) return(layered())
  if(nc == 1) return(layered(X))
  Y <- lapply(seq_len(nc), function(j,x) x[,j], x=X)
  names(Y) <- colnames(X)
  return(layered(LayerList=Y))
}

scalardilate.msr <- function(X, f, ...) {
  X$loc <- Xloc <- scalardilate(X$loc, f, ...)
  putlastshift(X, getlastshift(Xloc))
}
