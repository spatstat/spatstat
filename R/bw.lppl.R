#'
#'   bw.lppl.R
#' 
#'   Likelihood cross-validation for kernel smoother of point pattern on network
#'
#'   $Revision: 1.1 $ $Date: 2020/04/08 06:41:46 $
#'

bw.lppl <- function(X, ..., srange=NULL, ns=16, sigma=NULL,
                    weights=NULL, distance=c("euclidean", "path"),
                    shortcut=TRUE, warn=TRUE) {
  stopifnot(is.lpp(X))
  distance <- match.arg(distance)
  if(!is.null(sigma)) {
    stopifnot(is.numeric(sigma) && is.vector(sigma))
    ns <- length(sigma)
  } else {
    if(!is.null(srange)) check.range(srange) else {
      dd <- diameter(Frame(X))                                               
      ss <- bw.scott.iso(X)
      srange <- range(c(ss/10, ss*5, dd/5))
    }
    sigma <- sqrt(seq(from=srange[1L]^2, to=srange[2L]^2, length.out=ns))
  }
  a <- switch(distance,
              path = lcvlppHeat(X, sigma, ...,
                                weights=weights, shortcut=shortcut,
                                Nsave=max(128, 2*ns)),
              euclidean = lcvlppQuick(X, sigma, ...,
                                      weights=weights, shortcut=shortcut))
  result <- with(a,
                 bw.optim(cv, sigma, optimum="max",
                          creator="bw.lppl",
                          criterion="Likelihood Cross-Validation",
                          warnextreme=warn,
                          hargnames="srange",
                          unitname=unitname(X)))
  return(result)
}

lcvlppQuick <- function(X, sigmas, ..., weights=NULL, shortcut=TRUE) {
  ns <- length(sigmas)
  cv <- numeric(ns)
  if(shortcut) {
    #' omit calculation of integral term
    #' precompute the geometry data
    lam1 <- densityQuick.lpp(X, sigma=sigmas[1], weights=weights, ...,
                             savecomputed=TRUE)
    precooked <- attr(lam1, "savedstuff")
    for(i in 1:ns) {
      si <- sigmas[i]
      lamx <- densityQuick.lpp(X, sigma=si,
                               at="points", leaveoneout=TRUE,
                               weights=weights, 
                               precomputed=precooked, ...)
      lamx <- pmax(0, lamx)
      cv[i] <- sum(log(lamx))
    }
  } else {
    #' full calculation
    precooked <- NULL
    cooking <- TRUE
    for(i in 1:ns) {
      si <- sigmas[i]
      lamx <- densityQuick.lpp(X, sigma=si, at="points", leaveoneout=TRUE,
                               weights=weights, 
                               precomputed=precooked, ...)
      lam <- densityQuick.lpp(X, sigma=si,
                              weights=weights, 
                              precomputed=precooked,
                              savecomputed=cooking,
                              ...)
      if(cooking) {
        #' save geometry info for re-use in subsequent iterations
        precooked <- attr(lam, "savedstuff")
        cooking <- FALSE
      }
      lamx <- pmax(0, lamx)
      cv[i] <- sum(log(lamx)) - integral(lam)
    }
  }
  return(list(cv=cv, sigma=sigmas))
}

lcvlppHeat <- function(X, sigmas, ..., weights=NULL, shortcut=TRUE,
                       finespacing=FALSE, verbose=FALSE, Nsave=128,
                       dt=NULL,dx=NULL,iterMax=1e6) {
  #' first determine a resolution that is independent of sigmas
  #' for consistency between different calls.
  smax <- mean(sidelengths(Frame(X)))/4
  p <- resolve.heat.steps(smax, L=domain(X), dx=dx, dt=dt, iterMax=iterMax,
                          verbose=verbose)
  #' density at data points
  lamX <- densitypointsLPP(X, sigma=max(sigmas), ...,
                           dt=p$dt, dx=p$dx, finespacing=finespacing,
                           leaveoneout=TRUE, 
                           weights=weights, nsigma=Nsave,
                           verbose=verbose)
  #' map desired sigmas to actual calculated sigmas
  actualsigmas <- attr(lamX, "sigma")
  kleft <- findInterval(sigmas, actualsigmas,
                        all.inside=TRUE, rightmost.closed=TRUE)
  kright <- pmin(kleft+1L, Nsave)
  useleft <- (sigmas - actualsigmas[kleft]) < (actualsigmas[kright] - sigmas)
  kmap <- ifelse(useleft, kleft, kright)
  #' extract values for desired sigmas
  lamX <- lamX[ , kmap, drop=FALSE]
  lamX[] <- pmax(0, lamX)
  #' compute cross-validation term
  cv <- colSums(log(lamX))
  #'
  if(shortcut) {
    cv <- cv - npoints(X)
  } else {
    #' add integral term
    L <- as.linnet(X)
    lam <- densityfun(X, sigma=max(sigmas), ..., finespacing=finespacing,
                      weights=weights, nsigma=Nsave, verbose=verbose)
    for(j in seq_along(cv)) {
      kj <- kmap[j]
      lamjfun <- function(x, y, seg, tp) { lam(x, y, seg, tp, k=kj) }
      lamjim <- as.linim(linfun(lamjfun, L))
      cv[j] <- cv[j] - integral(lamjim)
    }
  }
  return(list(cv=cv, sigma=actualsigmas[kmap]))
}
