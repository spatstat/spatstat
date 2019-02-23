#
# Jinhom.R
#
#  $Revision: 1.11 $ $Date: 2017/06/05 10:31:58 $
#

Ginhom <- function(X, lambda=NULL, lmin=NULL,
                   ...,
                   sigma=NULL, varcov=NULL,
                   r=NULL, breaks=NULL,
                   ratio=FALSE, update = TRUE,
                   warn.bias=TRUE, savelambda=FALSE) {
  
  stopifnot(is.ppp(X))

  npts <- npoints(X)
  W <- as.owin(X)
  areaW <- area(W)
  miss.update <- missing(update)

  # determine 'r' values
  rmaxdefault <- rmax.rule("G", W, npts/areaW)
  breaks <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)
  if(!breaks$even)
    stop("r values must be evenly spaced")
  r <- breaks$r
  rmax <- breaks$max
  nr <- length(r)

  dangerous <- "lambda"
  danger <- TRUE
  
  # Intensity values at data points
  if(is.null(lambda)) {
    # No intensity data provided
    danger <- FALSE
    # Estimate density at points by leave-one-out kernel smoothing
    lamX <- density(X, ..., sigma=sigma, varcov=varcov,
                      at="points", leaveoneout=TRUE)
    lambdaX <- as.numeric(lamX)
    # negative or zero values are due to numerical error
    lambdaX <- pmax.int(lambdaX, .Machine$double.eps)
  } else {
    # lambda values provided
    if(is.im(lambda)) 
      lambdaX <- safelookup(lambda, X)
    else if(is.ppm(lambda) || is.kppm(lambda) || is.dppm(lambda)) {
          model <- lambda
          if(!update) {
            ## just use intensity of fitted model
            lambdaX <- predict(lambda, locations=X, type="trend")
          } else {
            ## re-fit model to data X
            model <-
              if(is.ppm(model)) update(model, Q=X) else update(model, X=X)
            lambdaX <- fitted(model, dataonly=TRUE)
            danger <- FALSE
            if(miss.update) 
              warn.once(key="Ginhom.update",
                        "The behaviour of Ginhom when lambda is a ppm object",
                        "has changed (in spatstat 1.37-0 and later).",
                        "See help(Ginhom)")
          }
        } else if(is.function(lambda)) 
      lambdaX <- lambda(X$x, X$y)
    else if(is.numeric(lambda) && is.vector(as.numeric(lambda))) {
      lambdaX <- lambda
      check.nvector(lambdaX, npts)
    } else stop(paste(sQuote("lambda"),
                      "should be a vector, a pixel image, or a function"))
    # negative values are illegal
    minX <- min(lambdaX)
    if(minX < 0)
      stop("Negative values of lambda were encountered at data points")
    if(minX == 0)
      stop("Zero values of lambda were encountered at data points")
  }
  # Minimum intensity
  if(!is.null(lmin)) {
    check.1.real(lmin)
    stopifnot(lmin >= 0)
    if(lmin >= min(lambdaX))
      stop("lmin must be smaller than all values of lambda")
  } else {
    # Compute minimum value over window
    if(is.null(lambda)) {
      # extract previously selected smoothing bandwidth
      sigma <- attr(lamX, "sigma")
      varcov <- attr(lamX, "varcov")
      # estimate density on a pixel grid and minimise
      lam <- density(X, ..., sigma=sigma, varcov=varcov, at="pixels")
      lmin <- min(lam)
      # negative or zero values may occur due to numerical error
      lmin <- max(lmin, .Machine$double.eps)
    } else {
      if(is.im(lambda)) 
        lmin <- min(lambda)
      else if(is.ppm(lambda) || is.kppm(lambda) || is.dppm(lambda)) 
        lmin <- min(predict(lambda))
      else if(is.function(lambda)) 
        lmin <- min(as.im(lambda, W))
      else if(is.numeric(lambda) && is.vector(as.numeric(lambda))) 
        lmin <- min(lambdaX)
    }
    if(lmin < 0)
      stop("Negative values of intensity encountered")
    # ensure lmin < lambdaX
    lmin <- min(lmin, 0.95 * min(lambdaX))
  }
  ## Compute intensity factor
  lratio <- lmin/lambdaX
  vv <- 1 - lratio
  if(warn.bias) {
    ra <- range(lratio)
    if(ra[1] < 1e-6 || ra[2] > 1 - 1e-6)
      warning(paste("Possible bias: range of values of lmin/lambdaX is",
                    prange(signif(ra, 5))),
              call.=FALSE)
  }
  ## sort data points in order of increasing x coordinate
  xx <- X$x
  yy <- X$y
  oX <- fave.order(xx)
  xord <- xx[oX]
  yord <- yy[oX]
  vord <- vv[oX]
  # compute local cumulative products
  z <- .C("locprod",
          n = as.integer(npts),
          x = as.double(xord),
          y = as.double(yord),
          v = as.double(vord),
          nr = as.integer(nr),
          rmax = as.double(rmax),
          ans = as.double(numeric(npts * nr)),
          PACKAGE = "spatstat")
  ans <- matrix(z$ans, nrow=nr, ncol=npts)
  # revert to original ordering
  loccumprod <- matrix(,  nrow=nr, ncol=npts)
  loccumprod[, oX] <- ans
  # border correction
  bX <- bdist.points(X)
  ok <- outer(r, bX, "<=")
  denom <- .rowSums(ok, nr, npts)
  loccumprod[!ok] <- 0
  numer <- .rowSums(loccumprod, nr, npts)
  # pack up
  Gdf <- data.frame(r=r, theo = 1 - exp(- lmin * pi * r^2))
  desc <- c("distance argument r", "theoretical Poisson %s")
  theo.denom <- rep.int(npts, nr)
  G <- ratfv(Gdf, NULL, theo.denom,
             "r", quote(G[inhom](r)),
             "theo", NULL, c(0,rmax),
             c("r", "{%s[%s]^{pois}}(r)"),
             desc,
             fname=c("G", "inhom"),
             ratio=ratio)
  G <- bind.ratfv(G,
                  data.frame(bord=denom-numer), denom,
                   "{hat(%s)[%s]^{bord}}(r)",
                  "border estimate of %s",
                  "bord",
                  ratio=ratio)
  # 
  formula(G) <- . ~ r
  fvnames(G, ".") <- c("bord", "theo")
  unitname(G) <- unitname(X)
  if(ratio)
    G <- conform.ratfv(G)
  if(danger)
    attr(G, "dangerous") <- dangerous
  if(savelambda) {
    attr(G, "lambda") <- lambdaX
    attr(G, "lmin") <- lmin
  }
  return(G)
}


   

Finhom <- function(X, lambda=NULL, lmin=NULL,
                   ...,
                   sigma=NULL, varcov=NULL,
                   r=NULL, breaks=NULL,
                   ratio=FALSE, update = TRUE,
                   warn.bias=TRUE, savelambda=FALSE) {
  
  stopifnot(is.ppp(X))

  npts <- npoints(X)
  W <- as.owin(X)
  areaW <- area(W)
  miss.update <- missing(update)

  # determine 'r' values
  rmaxdefault <- rmax.rule("F", W, npts/areaW)
  breaks <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)
  if(!breaks$even)
    stop("r values must be evenly spaced")
  r <- breaks$r
  rmax <- breaks$max
  nr <- length(r)

  dangerous <- "lambda"
  danger <- TRUE
  
  # Intensity values at data points
  if(is.null(lambda)) {
    # No intensity data provided
    danger <- FALSE
    # Estimate density at points by leave-one-out kernel smoothing
    lamX <- density(X, ..., sigma=sigma, varcov=varcov,
                      at="points", leaveoneout=TRUE)
    lambdaX <- as.numeric(lamX)
    # negative or zero values are due to numerical error
    lambdaX <- pmax.int(lambdaX, .Machine$double.eps)
  } else {
    # lambda values provided
    if(is.im(lambda)) 
      lambdaX <- safelookup(lambda, X)
    else if(is.ppm(lambda) || is.kppm(lambda) || is.dppm(lambda)) {
          model <- lambda
          if(!update) {
            ## just use intensity of fitted model
            lambdaX <- predict(lambda, locations=X, type="trend")
          } else {
            ## re-fit model to data X
            model <-
              if(is.ppm(model)) update(model, Q=X) else update(model, X=X)
            lambdaX <- fitted(model, dataonly=TRUE)
            danger <- FALSE
            if(miss.update) 
              warn.once(key="Finhom.update",
                        "The behaviour of Finhom when lambda is a ppm object",
                        "has changed (in spatstat 1.37-0 and later).",
                        "See help(Finhom)")
          }
        } else if(is.function(lambda)) 
      lambdaX <- lambda(X$x, X$y)
    else if(is.numeric(lambda) && is.vector(as.numeric(lambda))) {
      lambdaX <- lambda
      check.nvector(lambdaX, npts)
    } else stop(paste(sQuote("lambda"),
                      "should be a vector, a pixel image, or a function"))
    # negative values are illegal
    minX <- min(lambdaX)
    if(minX < 0)
      stop("Negative values of lambda were encountered at data points")
    if(minX == 0)
      stop("Zero values of lambda were encountered at data points")
  }
  # Minimum intensity
  if(!is.null(lmin)) {
    check.1.real(lmin)
    stopifnot(lmin >= 0)
    if(lmin >= min(lambdaX))
      stop("lmin must be smaller than all values of lambda")
  } else {
    # Compute minimum value over window
    if(is.null(lambda)) {
      # extract previously selected smoothing bandwidth
      sigma <- attr(lamX, "sigma")
      varcov <- attr(lamX, "varcov")
      # estimate density on a pixel grid and minimise
      lam <- density(X, ..., sigma=sigma, varcov=varcov, at="pixels")
      lmin <- min(lam)
      # negative or zero values may occur due to numerical error
      lmin <- max(lmin, .Machine$double.eps)
    } else {
      if(is.im(lambda)) 
        lmin <- min(lambda)
      else if(is.ppm(lambda) || is.kppm(lambda) || is.dppm(lambda)) 
        lmin <- min(predict(lambda))
      else if(is.function(lambda)) 
        lmin <- min(as.im(lambda, W))
      else if(is.numeric(lambda) && is.vector(as.numeric(lambda))) 
        lmin <- min(lambdaX)
    }
    if(lmin < 0)
      stop("Negative values of intensity encountered")
    # ensure lmin < lambdaX
    lmin <- min(lmin, 0.95 * min(lambdaX))
  }
  # Compute intensity factor
  lratio <- lmin/lambdaX
  vv <- 1 - lratio
  if(warn.bias) {
    ra <- range(lratio)
    if(ra[1] < 1e-6 || ra[2] > 1 - 1e-6)
      warning(paste("Possible bias: range of values of lmin/lambdaX is",
                    prange(signif(ra, 5))),
              call.=FALSE)
  }
  ## sort data points in order of increasing x coordinate
  xx <- X$x
  yy <- X$y
  oX <- fave.order(xx)
  xord <- xx[oX]
  yord <- yy[oX]
  vord <- vv[oX]
  # determine pixel grid and compute distance to boundary
  M <- do.call.matched(as.mask, append(list(w=W), list(...)))
  bM <- bdist.pixels(M, style="matrix")
  bM <- as.vector(bM)
  # x, y coordinates of pixels are already sorted by increasing x
  xM <- as.vector(rasterx.mask(M))
  yM <- as.vector(rastery.mask(M))
  nM <- length(xM)
  # compute local cumulative products
  z <- .C("locxprod",
         ntest = as.integer(nM),
         xtest = as.double(xM),
         ytest = as.double(yM),
         ndata = as.integer(npts),
         xdata = as.double(xord),
         ydata = as.double(yord),
         vdata = as.double(vord),
         nr = as.integer(nr),
         rmax = as.double(rmax),
         ans = as.double(numeric(nM * nr)),
         PACKAGE = "spatstat")
  loccumprod <- matrix(z$ans, nrow=nr, ncol=nM)
  # border correction
  ok <- outer(r, bM, "<=")
  denom <- .rowSums(ok, nr, nM)
  loccumprod[!ok] <- 0
  numer <- .rowSums(loccumprod, nr, nM)
  # pack up
  Fdf <- data.frame(r=r, theo = 1 - exp(- lmin * pi * r^2))
  desc <- c("distance argument r", "theoretical Poisson %s")
  theo.denom <- rep.int(npts, nr)
  FX <- ratfv(Fdf, NULL, theo.denom,
              "r",
              quote(F[inhom](r)),
              "theo", NULL, c(0,rmax),
              c("r","{%s[%s]^{pois}}(r)"),
              desc,
              fname=c("F", "inhom"),
              ratio=ratio)
  FX <- bind.ratfv(FX,
                  data.frame(bord=denom-numer), denom,
                  "{hat(%s)[%s]^{bord}}(r)",
                  "border estimate of %s",
                  "bord",
                  ratio=ratio)
  # 
  formula(FX) <- . ~ r
  fvnames(FX, ".") <- c("bord", "theo")
  unitname(FX) <- unitname(X)
  if(ratio)
    FX <- conform.ratfv(FX)
  if(danger)
    attr(FX, "dangerous") <- dangerous
  if(savelambda) {
    attr(FX, "lambda") <- lambdaX
    attr(FX, "lmin") <- lmin
  }
  return(FX)
}

Jinhom <- function(X, lambda=NULL, lmin=NULL,
                   ...,
                   sigma=NULL, varcov=NULL,
                   r=NULL, breaks=NULL, update = TRUE,
                   warn.bias=TRUE, savelambda=FALSE) {
  if(missing(update) & (is.ppm(lambda) || is.kppm(lambda) || is.dppm(lambda)))
    warn.once(key="Jinhom.update",
              "The behaviour of Jinhom when lambda is a ppm object",
              "has changed (in spatstat 1.37-0 and later).",
              "See help(Jinhom)")

  ## compute inhomogeneous G (including determination of r and lmin)
  GX <- Ginhom(X, lambda=lambda, lmin=lmin, ...,
               sigma=sigma, varcov=varcov, r=r, breaks=breaks,
               ratio=FALSE, update=update,
               warn.bias=warn.bias,
               savelambda=TRUE)
  ## extract auxiliary values to be used for Finhom
  r <- GX$r
  lmin <- attr(GX, "lmin")
  lambdaX <- attr(GX, "lambda")
  ## compute inhomogeneous J using previously-determined values
  FX <- Finhom(X, lambda=lambdaX, lmin=lmin, ...,
               sigma=sigma, varcov=varcov, r=r, ratio=FALSE, update=update,
               warn.bias=FALSE, savelambda=FALSE)
  ## evaluate inhomogeneous J function
  JX <- eval.fv((1-GX)/(1-FX))
  # relabel the fv object
  JX <- rebadge.fv(JX, quote(J[inhom](r)), c("J","inhom"),
                  names(JX), new.labl=attr(GX, "labl"))
  # tack on extra info
  attr(JX, "G") <- GX
  attr(JX, "F") <- FX
  attr(JX, "dangerous") <- attr(GX, "dangerous")
  if(savelambda) {
    attr(JX, "lmin") <- lmin
    attr(JX, "lambda") <- lambdaX
  }
  return(JX)
}
