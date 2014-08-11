#
#  leverage.R
#
#  leverage and influence
#
#  $Revision: 1.21 $  $Date: 2011/12/04 10:54:49 $
#

leverage <- function(model, ...) {
  UseMethod("leverage")
}

leverage.ppm <- function(model, ...,
                         drop=FALSE, iScore=NULL, iHessian=NULL, iArgs=list())
{
  fitname <- short.deparse(substitute(model))
  u <- list(fit=model, fitname=fitname)
  s <- ppm.influence(model, what="leverage", drop=drop,
                         iScore=iScore, iHessian=iHessian, iArgs=iArgs,
                         ...)
  a <- append(u, s)
  class(a) <- "leverage.ppm"
  return(a)
}

influence.ppm <- function(model, ...,
                          drop=FALSE, iScore=NULL, iHessian=NULL, iArgs=list())
{
  fitname <- short.deparse(substitute(model))
  u <- list(fit=model,fitname=fitname)
  s <- ppm.influence(model, what="influence", drop=drop,
                         iScore=iScore, iHessian=iHessian, iArgs=iArgs,
                         ...)
  a <- append(u, s)
  class(a) <- "influence.ppm"
  return(a)
}

dfbetas.ppm <- function(model, ..., drop=FALSE, iScore=NULL, iHessian=NULL, iArgs=list()) {
  fitname <- short.deparse(substitute(model))
  u <- list(fit=model,fitname=fitname)
  s <- ppm.influence(model, what="dfbetas", drop=drop,
                         iScore=iScore, iHessian=iHessian, iArgs=iArgs,
                         ...)
  a <- s$dfbetas
  attr(a, "info") <- u
  return(a)
}

ppm.influence <- function(fit,
           what=c("leverage", "influence", "dfbetas", "derivatives"),
           ..., iScore=NULL, iHessian=NULL, iArgs=list(), drop=FALSE,
           method=c("C", "interpreted")) {
  stopifnot(is.ppm(fit))
#  stopifnot(is.poisson.ppm(fit))
  what <- match.arg(what)
  method <- match.arg(method)
  #
  gotScore <- !is.null(iScore)
  gotHess <- !is.null(iHessian)
  needHess <- gotScore && any(what %in% c("leverage", "influence", "dfbetas"))
  if(!gotHess && needHess)
    stop("Must supply iHessian")
  #
  p <- length(coef(fit))
  vc <- vcov(fit, hessian=TRUE)
  fush <- hess <- solve(vc)
  lam <- fitted(fit, check=FALSE)
  mom <- model.matrix(fit)
  Q <- quad.ppm(fit)
  # hess = negative hessian of log (pseudo) likelihood
  # fush = E(hess)
  # invhess = solve(hess)
  # vc = solve(fush)
  #
  w <- w.quad(Q)
  loc <- union.quad(Q)
  isdata <- is.data(Q)
  #
  if(length(w) != length(lam))
    stop(paste("Internal error: length(w) = ", length(w),
               "!=", length(lam), "= length(lam)\n"))
  #
  if(!is.null(iScore)) {
    # evaluate additional (`irregular') components of score
    iscoredf <- mpl.get.covariates(iScore, loc, covfunargs=iArgs)
    iscoremat <- as.matrix(iscoredf)
    # count regular and irregular parameters
    nreg <- ncol(mom)
    nirr <- ncol(iscoremat)
    # add extra columns to model matrix
    mom <- cbind(mom, iscoremat)
    # evaluate additional (`irregular') entries of Hessian
    if(gotHess) {
      ihessdf <- mpl.get.covariates(iHessian, loc, covfunargs=iArgs)
      ihessmat <- as.matrix(ihessdf)
    }
    # recompute negative Hessian of log PL and its mean
    fush <- hessextra <- matrix(0, ncol(mom), ncol(mom))
    sub <- nreg + 1:nirr
    # integral over domain
    switch(method,
           interpreted = {
             for(i in seq(loc$n)) {
               # weight for integrand
               wti <- lam[i] * w[i]
               if(all(is.finite(wti))) {
                 # integral of outer product of score 
                 momi <- mom[i, ]
                 v1 <- outer(momi, momi, "*") * wti
                 if(all(is.finite(v1)))
                   fush <- fush + v1
                 # integral of Hessian
                 # contributions nonzero for irregular parameters
                 if(gotHess) {
                   v2 <- matrix(as.numeric(ihessmat[i,]), nirr, nirr) * wti
                   if(all(is.finite(v2)))
                     hessextra[sub, sub] <- hessextra[sub, sub] + v2
                 }
               }
             }
             # subtract sum over data points
             if(gotHess) {
               for(i in which(isdata)) {
                 v2 <- matrix(as.numeric(ihessmat[i,]), nirr, nirr) 
                 if(all(is.finite(v2)))
                   hessextra[sub, sub] <- hessextra[sub, sub] - v2
               }
               hess <- fush + hessextra
               invhess <- solve(hess)
             } else {
               invhess <- hess <- NULL
             }
           },
           C = {
             wlam <- lam * w
             fush <- sumouter(mom, wlam)
             if(gotHess) {
               # integral term
               ok <- is.finite(wlam) & apply(is.finite(ihessmat), 1, all)
               vintegral <-
                 if(all(ok)) wlam %*% ihessmat else
                             wlam[ok] %*% ihessmat[ok,, drop=FALSE]
               # sum over data points
               vdata <- colSums(ihessmat[isdata, , drop=FALSE], na.rm=TRUE)
               vcontrib <- vintegral - vdata
               hessextra[sub, sub] <-
                 hessextra[sub, sub] + matrix(vcontrib, nirr, nirr)
               hess <- fush + hessextra
               invhess <- solve(hess)
             } else {
               invhess <- hess <- NULL
             }
           })
    vc <- solve(fush)
  }
  if(!needHess) {
    hess <- fush
    invhess <- vc
  }
  #
  if(drop) {
    ok <- complete.cases(mom)
    Q <- Q[ok]
    mom <- mom[ok,]
    loc <- loc[ok]
    lam <- lam[ok]
    w   <- w[ok]
    isdata <- isdata[ok]
  }
  #
  result <- list()
  if("derivatives" %in% what) {
    rawresid <- isdata - lam * w
    score <- matrix(rawresid, nrow=1) %*% mom
    result$deriv <- list(mom=mom, score=score,
                         fush=fush, vc=vc,
                         hess=hess, invhess=invhess)
  }
  nloc <- npoints(loc)
  switch(method,
         interpreted = {
           b <- numeric(nloc)
           for(i in seq(nloc)) {
             momi <- mom[i,, drop=FALSE]
             b[i] <- momi %*% invhess %*% t(momi)
           }
         },
         C = {
           b <- quadform(mom, invhess)
         })
  if("leverage" %in% what) {
    # values of leverage (diagonal) at points of 'loc'
    h <- b * lam
    levval <- loc %mark% h
    levsmo <- smooth.ppp(levval, sigma=max(nndist(loc)))
    # nominal mean level
    a <- area.owin(loc$window)
    levmean <- p/a
    lev <- list(val=levval, smo=levsmo, ave=levmean)
    result$lev <- lev
  }
  if("influence" %in% what) {
    # values of influence at data points
    X <- loc[isdata]
    M <- (1/p) * b[isdata]
    V <- X %mark% M
    result$infl <- V
  }
  if("dfbetas" %in% what) {
    vex <- invhess %*% t(mom)
    switch(method,
           interpreted = {
             dis <- con <- matrix(0, nloc, ncol(mom))
             for(i in seq(nloc)) {
               vexi <- vex[,i, drop=FALSE]
               dis[i, ] <- isdata[i] * vexi
               con[i, ] <- - lam[i] * vexi
             }
           },
           C = {
             tvex <- t(vex)
             dis <- isdata * tvex
             con <- - lam  * tvex
           })
    colnames(dis) <- colnames(con) <- colnames(mom)
    result$dfbetas <- msr(Q, dis[isdata, ], con)
  }
  return(result)
}

plot.leverage.ppm <- function(x, ..., showcut=TRUE) {
  fitname <- x$fitname
  defaultmain <- paste("Leverage for", fitname)
  y <- x$lev
  do.call("plot.im",
          resolve.defaults(list(y$smo),
                           list(...),
                           list(main=defaultmain)))
  if(showcut) 
    contour(y$smo, levels=y$ave, add=TRUE, drawlabels=FALSE)
  invisible(NULL)
}

plot.influence.ppm <- function(x, ...) {
  fitname <- x$fitname
  defaultmain <- paste("Influence for", fitname)
  do.call("plot.ppp",
          resolve.defaults(list(x$infl),
                           list(...),
                           list(main=defaultmain)))
}

as.im.leverage.ppm <- function(X, ...) {
  return(X$lev$smo)
}

as.ppp.influence.ppm <- function(X, ...) {
  return(X$infl)
}
