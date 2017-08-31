#
#   measures.R
#
#  signed/vector valued measures with atomic and diffuse components
#
#  $Revision: 1.71 $  $Date: 2017/08/31 08:27:12 $
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
                atommass   = marksubset(data$discrete, data$atoms))
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
                     multiplot=TRUE,
                     massthresh=0,
                     equal.markscale=FALSE,
                     equal.ribbon=FALSE) {
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
  #' ensure image of density is present
  y <- solapply(y, augment.msr)

  #' ready to plot
  if(length(y) > 1) {
    ## plot as an array of panels
    userarg <- list(...)
    rowcol <- list(nrows=k, ncols=d)
    if(any(c("nrows", "ncols") %in% names(userarg))) rowcol <- list()
    #' determine common scales if required
    scaleinfo <- list()
    if(equal.markscale) {
      W <- Window(x)
      #' extract vectors of atomic masses from each panel
      marx <- lapply(y, with, "atommass")
      #' make a separate scale calculation for each panel
      scales <- sapply(marx, mark.scale.default, w=W, ...)
      scaleinfo$markscale <- min(scales)
      scaleinfo$markrange <- range(unlist(marx))
    } 
    if(equal.ribbon) {
      images <- lapply(y, attr, which="smoothdensity")
      scaleinfo$zlim <- range(sapply(images, range))
    } 
    ## go
    result <- do.call(plot.solist,
                      resolve.defaults(list(y),
                                       userarg,
                                       rowcol,
                                       scaleinfo,
                                       list(how=how,
                                            main=main,
                                            equal.scales=TRUE,
                                            halign=TRUE,
                                            valign=TRUE,
                                            claim.title.space=TRUE)))
    return(invisible(result))
  }
  ## scalar measure
  x <- y[[1]]
  ## get atoms
  xatomic <- (x$loc %mark% x$discrete)[x$atoms]
  if(length(massthresh) && all(is.finite(massthresh))) {
    ## ignore atoms with absolute mass <= massthresh
    check.1.real(massthresh)
    xatomic <- xatomic[abs(marks(xatomic)) > massthresh]
  }
  xtra.im <- graphicsPars("image")
  xtra.pp <- setdiff(graphicsPars("ppp"), c("box", "col"))
  xtra.pp <- union(xtra.pp, c("markrange", "marklevels"))
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
                        "drawlabels", "method", "vfont", "lty", "lwd",
                        "claim.title.space"))
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

Ops.msr <- function(e1,e2=NULL){
  vn <- c("val", "discrete", "density")
  if(nargs() == 1L) {
    #' unary operator
    if(!is.element(.Generic, c("+", "-")))
      stop(paste("Unary operation",
                 sQuote(paste0(.Generic, "A")),
                 "is undefined for a measure A."),
           call.=FALSE)
    e1 <- unclass(e1)
    e1[vn] <- lapply(e1[vn], .Generic)
    class(e1) <- "msr"
    return(e1)
  } else {
    #' binary operator
    m1 <- inherits(e1, "msr")
    m2 <- inherits(e2, "msr")
    if(m1 && m2) {
      if(!is.element(.Generic, c("+", "-")))
        stop(paste("Operation", sQuote(paste0("A", .Generic, "B")),
                   "is undefined for measures A, B"),
             call.=FALSE)
      k1 <- dim(e1)[2]
      k2 <- dim(e2)[2]
      if(k1 != k2) 
        stop(paste("Operation", sQuote(paste0("A", .Generic, "B")),
                   "is undefined because A, B have incompatible dimensions:",
                   "A is", ngettext(k1, "scalar", paste0(k1, "-vector")),
                   ", B is", ngettext(k2, "scalar", paste0(k2, "-vector"))),
             call.=FALSE)
      if(!identical(e1$loc, e2$loc)) {
        haha <- harmonise(e1, e2)
        e1 <- haha[[1L]]
        e2 <- haha[[2L]]
      }
      e1 <- unclass(e1)
      e2 <- unclass(e2)
      e1[vn] <- mapply(.Generic, e1[vn], e2[vn],
                       SIMPLIFY=FALSE)
      class(e1) <- "msr"
      return(e1)
    } else if(m1 && is.numeric(e2)) {
      if(!is.element(.Generic, c("/", "*")))
        stop(paste("Operation",
                   sQuote(paste0("A", .Generic, "z")),
                   "is undefined for a measure A and numeric z."),
             call.=FALSE)
      e1 <- unclass(e1)
      e1[vn] <- lapply(e1[vn], .Generic, e2=e2)
      class(e1) <- "msr"
      return(e1)
    } else if(m2 && is.numeric(e1)) {
      if(.Generic != "*") 
        stop(paste("Operation",
                   sQuote(paste0("z", .Generic, "A")),
                   "is undefined for a measure A and numeric z."),
             call.=FALSE)
      e2 <- unclass(e2)
      e2[vn] <- lapply(e2[vn], .Generic, e1=e1)
      class(e2) <- "msr"
      return(e2)
    }
    stop(paste("Operation", sQuote(paste0("e1", .Generic, "e2")),
               "is undefined for this kind of data"),
         call.=FALSE)
  }
}

measurePositive <- function(x) {
  if(!inherits(x, "msr"))
    stop("x must be a measure", call.=FALSE)
  y <- x
  y$discrete <- pmax(0, x$discrete)
  y$density  <- pmax(0, x$density)
  y$val      <- y$discrete + y$wt * y$density
  return(y)
}

measureNegative <- function(x) {
  if(!inherits(x, "msr"))
    stop("x must be a measure", call.=FALSE)
  y <- x
  y$discrete <- -pmin(0, x$discrete)
  y$density  <- -pmin(0, x$density)
  y$val      <- y$discrete + y$wt * y$density
  return(y)
}

measureVariation <- function(x) {
  if(!inherits(x, "msr"))
    stop("x must be a measure", call.=FALSE)
  y <- x
  y$discrete <- abs(x$discrete)
  y$density  <- abs(x$density)
  y$val      <- y$discrete + y$wt * y$density
  return(y)
}

totalVariation <- function(x) integral(measureVariation(x))
  
harmonise.msr <- local({

  harmonise.msr <- function(...) {
    argz <- list(...)
    n <- length(argz)
    if(n == 0) return(argz)
    ismeasure <- sapply(argz, inherits, what="msr")
    if(!any(ismeasure))
      stop("No measures supplied")
    if(!all(ismeasure))
    stop("All arguments should be measures (objects of class msr)")
    if(n < 2) return(argz)
    result <- vector(mode="list", length=n)
    ## extract entries
    loclist <- lapply(argz, getElement, name="loc")
    atomlist <- lapply(argz, getElement, name="atoms")
    masslist <- lapply(argz, getElement, name="discrete")
    denslist <- lapply(argz, getElement, name="density")
    ## check for compatible dimensions of measure values
    dimen <- unique(sapply(argz, ncol))
    if(length(dimen) > 1)
      stop("Measures have different dimensions:", commasep(sort(dimen)))
    ## check for marked points
    ismarked <- sapply(loclist, is.marked)
    if(any(ismarked) && !all(ismarked))
      stop("Some, but not all, quadrature schemes are marked")
    ismarked <- all(ismarked)
    ## union of all quadrature points in all measures
    Uloc <- do.call(superimpose, append(unname(loclist), list(check=FALSE)))
    Uloc <- unique(Uloc)
    nU <- npoints(Uloc)
    ## match each quadrature set to the union
    ## and find nearest data point to each point in the union
    if(!ismarked) {
      matchlist <- lapply(loclist, nncross, Y=Uloc, what="which")
      nearlist  <- lapply(loclist, ssorcnn, xx=Uloc, what="which")
    } else {
      stop("Not yet implemented for marked quadrature schemes")
    }
    ## nearest neighbour interpolation of density values of each argument
    ## onto the common quadrature set
    Udenslist <- mapply(extract, x=denslist, i=nearlist,
                        SIMPLIFY=FALSE)
    ## initialise other bits
    noatoms  <- logical(nU) 
    zeromass <- if(dimen == 1) numeric(nU) else matrix(0, nU, dimen)
    Uatomlist <- rep(list(noatoms), n)  
    Umasslist <- rep(list(zeromass), n)
    ## assign atoms in each argument
    Uatomlist <- mapply(subsetgets, x=Uatomlist, i=matchlist, value=atomlist,
                        SIMPLIFY=FALSE)
    Umasslist <- mapply(subsetgets, x=Umasslist, i=matchlist, value=masslist,
                        SIMPLIFY=FALSE)
    ## union of atoms
    isatom <- Reduce("|", Uatomlist)
    ## masses at atoms
    Umasslist <- lapply(Umasslist, extract, i=isatom)
    ## make common quadrature scheme
    UQ <- quadscheme(Uloc[isatom], Uloc[!isatom])
    ## reorder density data correspondingly
    neworder <- c(which(isatom), which(!isatom))
    Udenslist <- lapply(Udenslist, extract, i=neworder)
    ## make new measures
    result <- mapply(msr,
                     MoreArgs=list(qscheme=UQ),
                     discrete=Umasslist,
                     density=Udenslist,
                     SIMPLIFY=FALSE)
    names(result) <- names(argz)
    class(result) <- unique(c("solist", class(result)))
    return(result)
  }

  ssorcnn <- function(xx, yy, what) nncross(xx, yy, what=what)
  
  extract <- function(x, i) {
    if(is.matrix(x)) x[i, , drop=FALSE] else x[i]
  }
  subsetgets <- function(x, i, value) {
    if(is.matrix(x)) {
      x[i, ] <- value
    } else {
      x[i] <- value
    }
    return(x)
  }

  harmonise.msr
})
