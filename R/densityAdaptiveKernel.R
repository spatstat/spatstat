#'
#'   densityAdaptiveKernel.R
#'
#'   $Revision: 1.4 $  $Date: 2019/03/07 02:58:08 $
#'
#'
#'  Adaptive kernel smoothing via 3D FFT
#'

densityAdaptiveKernel <- function(X, ...) {
  UseMethod("densityAdaptiveKernel")
}

densityAdaptiveKernel.ppp <- function(X, bw, ...,
                                      weights=NULL,
                                      at=c("pixels", "points"),
                                      edge=TRUE, 
                                      ngroups) {
  stopifnot(is.ppp(X))
  at <- match.arg(at)
  nX <- npoints(X)

  if(nX == 0)
    switch(at,
           points = return(numeric(nX)),
           pixels = return(as.im(0, W=Window(X), ...)))
                     
  if(missing(ngroups) || is.null(ngroups)) {
    ngroups <- max(1L, floor(sqrt(npoints(X))))
  } else if(any(is.infinite(ngroups))) {
    ngroups <- nX
  } else {
    check.1.integer(ngroups)
    ngroups <- min(nX, ngroups)
  }

  if(weighted <- !is.null(weights)) {
    check.nvector(weights, nX, oneok=TRUE)
    if(length(weights) == 1) weights <- rep(weights, nX)
  } else weights <- rep(1,nX)

  ## determine bandwidth for each data point
  if(missing(bw) || is.null(bw)) {
    bw <- do.call.matched(bw.abram,
                          resolve.defaults(list(X=X, at="points"),
                                           list(...)),
                          extrargs=names(args(as.mask)))
  } else if(is.numeric(bw)) {
    check.nvector(bw, nX, oneok=TRUE)
    if(length(bw) == 1) bw <- rep(bw, nX)
  } else if(is.im(bw)) {
    bw <- safelookup(bw, X, warn=FALSE)
    if(anyNA(bw))
      stop("Some data points lie outside the domain of image 'bw'",
           call.=FALSE)
  } else if(inherits(bw, "funxy")) {
    bw <- bw(X)
    if(anyNA(bw))
      stop("Some data points lie outside the domain of function 'bw'",
           call.=FALSE)
  } else stop("Argument 'bw' should be a numeric vector or a pixel image")

  #' divide bandwidths into groups
  p <- seq(0,1,length=ngroups+1)
  qbands <- quantile(bw, p)
  groupid <- findInterval(bw,qbands,all.inside=TRUE)
  #' map to middle of group
  pmid <- (p[-1] + p[-length(p)])/2
  qmid   <- quantile(bw, pmid)

  marks(X) <- if(weighted) weights else NULL
  group <- factor(groupid, levels=1:ngroups)
  Y <- split(X, group)

  Z <- mapply(density.ppp,
              x=Y,
              sigma=as.list(qmid),
              weights=lapply(Y, marks),
              MoreArgs=list(edge=edge, at=at, ...),
              SIMPLIFY=FALSE)

  ZZ <- switch(at,
               pixels = im.apply(Z, "sum"),
               points = unsplit(Z, group))
  return(ZZ)
}

