#
# Jinhom.R
#
#  $Revision: 1.9 $ $Date: 2014/10/24 00:22:30 $
#

Ginhom <- function(X, lambda=NULL, lmin=NULL,
                   ...,
                   sigma=NULL, varcov=NULL,
                   r=NULL, breaks=NULL,
                   ratio=FALSE, update = TRUE) {
  
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
    lmin <- min(lmin, lambdaX)
  }
  # Compute intensity factor
  lratio <- lmin/lambdaX
  vv <- 1 - lratio
  bad <- (lratio > 1)
  if((nbad <- sum(bad)) > 0)
    stop(paste("Value of", sQuote("lmin"), "exceeds",
               nbad, gettext(nbad, "value", "values"),
               "of", sQuote("lambda")))
   # sort data points in order of increasing x coordinate
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
          ans = as.double(numeric(npts * nr)))
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
  return(G)
}


   

Finhom <- function(X, lambda=NULL, lmin=NULL,
                   ...,
                   sigma=NULL, varcov=NULL,
                   r=NULL, breaks=NULL,
                   ratio=FALSE, update = TRUE) {
  
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
    lmin <- min(lmin, lambdaX)
  }
  # Compute intensity factor
  lratio <- lmin/lambdaX
  vv <- 1 - lratio
  bad <- (lratio > 1)
  if((nbad <- sum(bad)) > 0)
    stop(paste("Value of", sQuote("lmin"), "exceeds",
               nbad, gettext(nbad, "value", "values"),
               "of", sQuote("lambda")))
  # sort data points in order of increasing x coordinate
  xx <- X$x
  yy <- X$y
  oX <- fave.order(xx)
  xord <- xx[oX]
  yord <- yy[oX]
  vord <- vv[oX]
  # determine pixel grid and compute distance to boundary
  M <- do.call.matched("as.mask", append(list(w=W), list(...)))
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
         ans = as.double(numeric(nM * nr)))
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
  return(FX)
}

Jinhom <- function(X, lambda=NULL, lmin=NULL,
                   ...,
                   sigma=NULL, varcov=NULL,
                   r=NULL, breaks=NULL, update = TRUE) {
  if(missing(update) & (is.ppm(lambda) || is.kppm(lambda) || is.dppm(lambda)))
    warn.once(key="Jinhom.update",
              "The behaviour of Jinhom when lambda is a ppm object",
              "has changed (in spatstat 1.37-0 and later).",
              "See help(Jinhom)")
        
  GX <- Ginhom(X, lambda=lambda, lmin=lmin, ...,
               sigma=sigma, varcov=varcov, r=r, breaks=breaks, ratio=FALSE, update=update)
  r <- GX$r
  FX <- Finhom(X, lambda=lambda, lmin=lmin, ...,
               sigma=sigma, varcov=varcov, r=r, ratio=FALSE, update=update)
  JX <- eval.fv((1-GX)/(1-FX))
  # relabel the fv object
  JX <- rebadge.fv(JX, quote(J[inhom](r)), c("J","inhom"),
                  names(JX), new.labl=attr(GX, "labl"))
  # tack on extra info
  attr(JX, "G") <- GX
  attr(JX, "F") <- FX
  attr(JX, "dangerous") <- attr(GX, "dangerous")
  return(JX)
}
