#
#   measures.R
#
#  signed/vector valued measures with atomic and diffuse components
#
#  $Revision: 1.60 $  $Date: 2016/06/30 03:57:09 $
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
  stopifnot(is.numeric(discrete) || is.logical(discrete))
  stopifnot(is.numeric(density))
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
  stuff <- list(increment  = data$val,
                is.atom    = data$atoms,
                discrete   = data$discrete,
                density    = data$density,
                continuous = data$density * data$wt,
                qweights   = data$wt,
                qlocations = data$loc,
                atoms      = data$loc[data$atoms],
                atommass   = data$wt[data$atoms])
  y <- eval(substitute(expr), envir=stuff, enclos=parent.frame())
  if(is.character(y) && length(y) == 1 && y %in% names(stuff))
    y <- stuff[[y]]
  return(y)
}

print.msr <- function(x, ...) {
  xloc <- x$loc
  n <- npoints(xloc)
  d <- ncol(as.matrix(x$val))
  splat(paste0(if(d == 1) "Scalar" else paste0(d, "-dimensional vector"),
               "-valued measure"))
  if(d > 1 && !is.null(cn <- colnames(x$val)) && waxlyrical("space"))
    splat("vector components:", commasep(sQuote(cn)))
  if(is.marked(xloc)) {
    splat("\tDefined on 2-dimensional space x marks")
    if(is.multitype(xloc))
      exhibitStringList("\tPossible marks: ", levels(marks(xloc)))
  } 
  if(waxlyrical("gory")) {
    splat("Approximated by", n, "quadrature points")
    print(as.owin(xloc))
    splat(sum(x$atoms), "atoms")
  }
  if(waxlyrical("extras")) {
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

is.multitype.msr <- function(X, ...) {
  is.multitype(X$loc, ...)
}
is.marked.msr <- function(X, ...) {
  is.marked(X$loc, ...)
}

split.msr <- function(x, f, drop=FALSE, ...) {
  xloc <- x$loc
  ## determine split using rules for split.ppp
  locsplit <- if(missing(f))
    split(xloc, drop=drop) else split(xloc, f, drop=drop)
  ## extract grouping factor 
  g <- attr(locsplit, "fgroup")
  ## split contributions to measure
  atomsplit <- split(x$atoms, g, drop=drop) # hyuk
  wtsplit <- split(x$wt, g, drop=drop)
  if(ncol(x) == 1) {
    ## scalar measure
    valsplit  <- split(x$val, g, drop=drop)
    discsplit <- split(x$discrete, g, drop=drop)
    denssplit <- split(x$density, g, drop=drop)
  } else {
    ## vector measure
    valsplit  <- lapply(split(as.data.frame(x$val), g, drop=drop),
                        as.matrix)
    discsplit <- lapply(split(as.data.frame(x$discrete), g, drop=drop),
                        as.matrix)
    denssplit <- lapply(split(as.data.frame(x$density), g, drop=drop),
                        as.matrix)
  }
  ## form the component measures
  result <- mapply(list,
                   loc=locsplit,
                   val=valsplit,
                   atoms=atomsplit,
                   discrete=discsplit,
                   density=denssplit,
                   wt=wtsplit,
                   SIMPLIFY=FALSE)
  names(result) <- names(locsplit)
  result <- lapply(result, "class<-", value="msr")
  if(drop && any(isnul <- (sapply(locsplit, npoints) == 0)))
    result[isnul] <- NULL
  result <- as.solist(result)
  return(result)
}

integral.msr <- function(f, domain=NULL, ...) {
  stopifnot(inherits(f, "msr"))
  if(!is.null(domain)) {
    if (is.tess(domain)) 
      return(sapply(tiles(domain), integral.msr, f = f))
    f <- f[domain]
  }
  y <- with(f, "increment")
  if(is.matrix(y)) apply(y, 2, sum) else sum(y)
}

augment.msr <- function(x, ..., sigma) {
  ## add a pixel image of the smoothed density component
  stopifnot(inherits(x, "msr"))
  if(!is.null(attr(x, "smoothdensity"))) return(x)
  d <- ncol(as.matrix(x$val))
  xloc <- x$loc
  W <- as.owin(xloc)
  if(missing(sigma)) sigma <- maxnndist(xloc, positive=TRUE)
  if(is.multitype(xloc)) {
    ## multitype case - split by type, smooth, sum
    y <- lapply(split(x), augment.msr, sigma=sigma, ...)
    z <- lapply(y, attr, which="smoothdensity")
    if((nc <- ncol(x)) == 1) {
      ## scalar valued
      smo <- Reduce("+", z)
    } else {
      ## vector valued
      smo <- vector(mode="list", length=nc)
      for(j in 1:nc) 
        smo[[j]] <- Reduce("+", lapply(z, "[[", i=j))
      smo <- as.solist(smo)
    }
    attr(x, "smoothdensity") <- smo
    return(x)
  }   
  ## smooth density unless constant
  xdensity <- as.matrix(x$density)
  ra <- apply(xdensity, 2, range)
  varble <- apply(as.matrix(ra), 2, diff) > sqrt(.Machine$double.eps)
  ##
  if(d == 1) {
    smo <- if(!varble) as.im(mean(xdensity), W=W) else
           do.call(Smooth,
                   resolve.defaults(list(X=xloc %mark% xdensity),
                                    list(...),
                                    list(sigma=sigma)))
  } else {
    smo <- vector(mode="list", length=d)
    names(smo) <- colnames(x)
    if(any(varble)) 
      smo[varble] <-
        do.call(Smooth,
                resolve.defaults(list(X=xloc %mark% xdensity[,varble, drop=FALSE]),
                                 list(...),
                                 list(sigma=sigma)))
    if(any(!varble)) 
      smo[!varble] <- lapply(apply(xdensity[, !varble, drop=FALSE], 2, mean),
                             as.im, W=W)
    smo <- as.solist(smo)
  }
  attr(x, "smoothdensity") <- smo
  return(x)
}

plot.msr <- function(x, ..., add=FALSE,
                     how=c("image", "contour", "imagecontour"),
                     main=NULL, 
                     do.plot=TRUE,
                     multiplot=TRUE) {
  if(is.null(main)) 
    main <- short.deparse(substitute(x))
  how <- match.arg(how)
  
  if(!multiplot) {
    ## compress everything to a single panel
    x$loc <- unmark(x$loc)
    if(is.matrix(x$val))      x$val <- rowSums(x$val)
    if(is.matrix(x$discrete)) x$discrete <- rowSums(x$discrete)
    if(is.matrix(x$density))  x$density <- rowSums(x$density)
    if(!is.null(smo <- attr(x, "smoothdensity")) && inherits(smo, "solist"))
      attr(x, "smoothdensity") <- Reduce("+", smo)
  }

  d <- dim(x)[2]
  k <- if(is.multitype(x)) length(levels(marks(x$loc))) else 1

  ## multiple plot panels may be generated
  if(k == 1 && d == 1) {
    ## single plot
    y <- solist(x)
  } else if(k > 1 && d == 1) {
    ## multitype
    y <- split(x)
  } else if(k == 1 && d > 1) {
    ## vector-valued
    y <- unstack(x)
  } else if(k > 1 && d > 1) {
    ## both multitype and vector-valued
    y <- split(x)
    typenames <- names(y)
    vecnames <- colnames(x$val)
    y <- as.solist(Reduce(append, lapply(y, unstack)))
    names(y) <- as.vector(t(outer(typenames, vecnames, paste, sep=".")))
  } 
  # ensure image of density is present
  y <- lapply(y, augment.msr)

  if(length(y) > 1) {
    ## plot as an array of panels
    result <- do.call(plot.solist, resolve.defaults(list(y),
                                                    list(...),
                                                    list(nrows=k, ncols=d),
                                                    list(how=how,
                                                         main=main,
                                                         equal.scales=TRUE)))
    return(invisible(result))
  }
  ## scalar measure
  x <- y[[1]]
  xatomic <- (x$loc %mark% x$discrete)[x$atoms]
  xtra.im <- graphicsPars("image")
  xtra.pp <- setdiff(graphicsPars("ppp"), c("box", "col"))
  xtra.ow <- graphicsPars("owin")
  smo <- attr(x, "smoothdensity")
  ##
  do.image <-  how %in% c("image", "imagecontour")
  do.contour <-  how %in% c("contour", "imagecontour")
  ## allocate space for plot and legend using do.plot=FALSE mechanism
  pdata <- do.call.matched(plot.ppp,
                           resolve.defaults(list(x=xatomic,
                                                 do.plot=FALSE,
                                                 main=main),
                                            list(...),
                                            list(show.all=TRUE)),
                           extrargs=xtra.pp)
  result <- pdata
  bb <- attr(pdata, "bbox")
  if(do.image) {
    idata <- do.call.matched(plot.im,
                             resolve.defaults(list(x=smo,
                                                   main=main,
                                                   do.plot=FALSE),
                                              list(...)),
                             extrargs=xtra.im)
    result <- idata
    bb <- boundingbox(bb, attr(idata, "bbox"))
  }
  ##
  attr(result, "bbox") <- bb
  ##
  if(do.plot) {
    if(!add) {
      blankmain <- prepareTitle(main)$blank
      ## initialise plot
      do.call.matched(plot.owin,
                      resolve.defaults(list(x=bb, type="n", main=blankmain),
                                       list(...)),
                      extrargs=xtra.ow)
    }
    ## display density
    if(do.image) 
      do.call.matched(plot.im,
                      resolve.defaults(list(x=smo, add=TRUE),
                                       list(...),
                                       list(main=main, show.all=TRUE)),
                      extrargs=xtra.im)
    if(do.contour) 
      do.call.matched(contour.im,
                      resolve.defaults(list(x=smo, add=TRUE),
                                       list(...),
                                       list(main=main,
                                            axes=FALSE, show.all=!do.image)),
                      extrargs=c("zlim", "labels", "labcex",
                        ## DO NOT ALLOW 'col' 
                        "drawlabels", "method", "vfont", "lty", "lwd"))
    ## display atoms
    do.call.matched(plot.ppp,
                    resolve.defaults(list(x=xatomic, add=TRUE, main=""),
                                     list(...),
                                     list(show.all=TRUE)),
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
    loci  <- locn[i, clip=TRUE]
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
    result <- solist(result)
  return(result)
}

as.owin.msr <- function(W, ..., fatal=TRUE) {
  as.owin(W$loc, ..., fatal=fatal)
}

domain.msr <- Window.msr <- function(X, ...) { as.owin(X) } 

shift.msr <- function(X,  ...) {
  X$loc <- Xloc <- shift(X$loc, ...)
  if(!is.null(smo <- attr(X, "smoothdensity")))
    attr(X, "smoothdensity") <- shift(smo, getlastshift(Xloc))
  putlastshift(X, getlastshift(Xloc))
}

as.layered.msr <- local({

  as.layered.msr <- function(X) {
    nc <- ncol(X)
    if(nc == 0) return(layered())
    if(nc == 1) return(layered(X))
    Y <- lapply(seq_len(nc), pickcol, x=X)
    names(Y) <- colnames(X)
    return(layered(LayerList=Y))
  }

  pickcol <- function(j,x) x[,j]
  
  as.layered.msr
})


scalardilate.msr <- function(X, f, ...) {
  X$loc <- Xloc <- scalardilate(X$loc, f, ...)
  putlastshift(X, getlastshift(Xloc))
}
