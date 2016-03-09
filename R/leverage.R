#
#  leverage.R
#
#  leverage and influence
#
#  $Revision: 1.54 $  $Date: 2016/03/09 10:31:14 $
#

leverage <- function(model, ...) {
  UseMethod("leverage")
}

leverage.ppm <- function(model, ...,
                         drop=FALSE, iScore=NULL, iHessian=NULL, iArgs=NULL)
{
  fitname <- short.deparse(substitute(model))
  u <- list(fitname=fitname, fit.is.poisson=is.poisson(model))
  s <- ppmInfluence(model, what="leverage", drop=drop,
                         iScore=iScore, iHessian=iHessian, iArgs=iArgs,
                         ...)
  a <- append(u, s)
  class(a) <- "leverage.ppm"
  return(a)
}

influence.ppm <- function(model, ...,
                          drop=FALSE, iScore=NULL, iHessian=NULL, iArgs=NULL)
{
  fitname <- short.deparse(substitute(model))
  u <- list(fitname=fitname, fit.is.poisson=is.poisson(model))
  s <- ppmInfluence(model, what="influence", drop=drop,
                         iScore=iScore, iHessian=iHessian, iArgs=iArgs,
                         ...)
  a <- append(u, s)
  class(a) <- "influence.ppm"
  return(a)
}

dfbetas.ppm <- function(model, ...,
                        drop=FALSE, iScore=NULL, iHessian=NULL, iArgs=NULL) {
  fitname <- short.deparse(substitute(model))
  u <- list(fitname=fitname, fit.is.poisson=is.poisson(model))
  s <- ppmInfluence(model, what="dfbetas", drop=drop,
                         iScore=iScore, iHessian=iHessian, iArgs=iArgs,
                     ...)
  a <- s$dfbetas
  attr(a, "info") <- u
  return(a)
}


ppmInfluence <- function(fit,
                         what=c("leverage", "influence", "dfbetas",
                                "derivatives", "increments"),
                         ...,
                         iScore=NULL, iHessian=NULL, iArgs=NULL,
                         drop=FALSE,
                         method=c("C", "interpreted"),
                          precomputed=list(),
                         sparseOK=FALSE) {
  stopifnot(is.ppm(fit))
  what <- match.arg(what, several.ok=TRUE)
  method <- match.arg(method)
  if(is.null(iArgs))
    iArgs <- fit$covfunargs
  gotScore <- !is.null(iScore)
  gotHess <- !is.null(iHessian)
  influencecalc <- any(what %in% c("leverage", "influence", "dfbetas"))
  needHess <- gotScore && influencecalc
  if(!gotHess && needHess)
    stop("Must supply iHessian")
  if(fit$method == "logi" && !spatstat.options("allow.logi.influence"))
    stop("ppm influence measures are not yet implemented for method=logi")
  sparse <- sparseOK && spatstat.options('developer')
  #
  # extract precomputed values if given
  theta  <- precomputed$coef   %orifnull% coef(fit)
  lam    <- precomputed$lambda %orifnull% fitted(fit, check=FALSE)
  mom    <- precomputed$mom    %orifnull% model.matrix(fit)
  # 
#  p <- length(theta)
  vc <- vcov(fit, hessian=TRUE)
  fush <- hess <- solve(vc)
  Q <- quad.ppm(fit)
  # hess = negative hessian of log (pseudo) likelihood
  # fush = E(hess)
  # invhess = solve(hess)
  # vc = solve(fush)
  #
  w <- w.quad(Q)
#  loc <- union.quad(Q)
#  isdata <- is.data(Q)
  #
  if(length(w) != length(lam))
    stop(paste("Internal error: length(w) = ", length(w),
               "!=", length(lam), "= length(lam)\n"))

  a <- ppmInfluenceEngine(fit      = fit,
                          what     = what,
                          iScore   = iScore,
                          iHessian = iHessian,
                          iArgs    = iArgs,
                          drop     = drop,
                          method   = method,
                          sparse   = sparse,
                          vc       = vc,
                          fush     = fush,
                          hess     = hess,
                          Q        = Q,
                          theta    = theta,
                          lam      = lam,
                          mom      = mom)
  return(a)
}

ppmInfluenceEngine <-
  function(fit, what, iScore, iHessian, iArgs,
           drop, method, sparse,
           vc, fush, hess,
           Q, theta, lam, mom) {
  gotScore <- !is.null(iScore)
  gotHess <- !is.null(iHessian)
  influencecalc <- any(what %in% c("leverage", "influence", "dfbetas"))
  needHess <- gotScore && influencecalc
  p <- length(theta)
  w <- w.quad(Q)
  loc <- union.quad(Q)
  isdata <- is.data(Q)
  ## evaluate additional (`irregular') components of score
  iscoremat <- ppmDerivatives(fit, "gradient", iScore, loc, covfunargs=iArgs)
  gotScore <- !is.null(iscoremat)
  needHess <- gotScore && influencecalc
  if(gotScore) {
    ## count regular and irregular parameters
    nreg <- ncol(mom)
    nirr <- ncol(iscoremat)
    ## add extra columns to model matrix
    mom <- cbind(mom, iscoremat)
    REG <- 1:nreg
    IRR <- nreg + 1:nirr
    ## evaluate additional (`irregular') entries of Hessian
    ihessmat <- ppmDerivatives(fit, "hessian", iHessian, loc, covfunargs=iArgs)
    if(gotHess <- !is.null(ihessmat))
    ## recompute negative Hessian of log PL and its mean
    fush <- hessextra <- matrix(0, ncol(mom), ncol(mom))
    ## integral over domain
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
                     hessextra[IRR, IRR] <- hessextra[IRR, IRR] + v2
                 }
               }
             }
             # subtract sum over data points
             if(gotHess) {
               for(i in which(isdata)) {
                 v2 <- matrix(as.numeric(ihessmat[i,]), nirr, nirr) 
                 if(all(is.finite(v2)))
                   hessextra[IRR, IRR] <- hessextra[IRR, IRR] - v2
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
               isfin <- is.finite(wlam) & apply(is.finite(ihessmat), 1, all)
               vintegral <-
                 if(all(isfin)) wlam %*% ihessmat else
                             wlam[isfin] %*% ihessmat[isfin,, drop=FALSE]
               # sum over data points
               vdata <- .colSums(ihessmat[isdata, , drop=FALSE],
                                 sum(isdata), ncol(ihessmat),
                                 na.rm=TRUE)
               vcontrib <- vintegral - vdata
               hessextra[IRR, IRR] <-
                 hessextra[IRR, IRR] + matrix(vcontrib, nirr, nirr)
               hess <- fush + hessextra
               invhess <- solve(hess)
             } else {
               invhess <- hess <- NULL
             }
           })
    vc <- solve(fush)
  } else {
    REG <- 1:ncol(mom)
  }
  
  if(!needHess) {
    hess <- fush
    invhess <- vc
  }
  #
  ok <- NULL
  if(drop) {
    ok <- complete.cases(mom)
    if(all(ok)) {
      ok <- NULL
    } else {
      Q <- Q[ok]
      mom <- mom[ok, , drop=FALSE]
      loc <- loc[ok]
      lam <- lam[ok]
      w   <- w[ok]
      isdata <- isdata[ok]
    }
  } 
  
  # second order interaction terms
  # ddS[i,j, ] = Delta_i Delta_j S(x)
  ddS <- NULL
  if(!all(what == "derivatives") && !is.poisson(fit)) {
    ddS <- deltasuffstat(fit, dataonly=FALSE, quadsub=ok, sparseOK=sparse)
    if(is.null(ddS)) {
      warning("Second order interaction terms are not implemented for this model; they are treated as zero")
    } else {
      sparse <- inherits(ddS, "sparse3Darray")
      if(gotScore) {
        ## add extra planes of zeroes to second-order model matrix
        ## (zero because the irregular components are part of the trend)
        paddim <- c(dim(ddS)[1:2], nirr)
        if(!sparse) {
          ddS <- abind::abind(ddS, array(0, dim=paddim), along=3)
        } else {
          ddS <- bind.sparse3Darray(ddS,
                                    sparse3Darray(dims=paddim),
                                    along=3)
        }
      }
    }
  }

  # ........  start assembling results .....................
  # 
  result <- list()
  # 
  if("derivatives" %in% what) {
    rawresid <- isdata - lam * w
    score <- matrix(rawresid, nrow=1) %*% mom
    result$deriv <- list(mom=mom, score=score,
                         fush=fush, vc=vc,
                         hess=hess, invhess=invhess)
  }
  if(all(what == "derivatives"))
    return(result)

  # compute effect of adding/deleting each quadrature point
  #    columns index the point being added/deleted
  #    rows index the points affected
  eff <- mom
  if(!is.poisson(fit) && !is.null(ddS)) {
    # effect of addition/deletion of U[j] on score contribution from data points
    ddSX <- ddS[isdata, , , drop=FALSE]
    eff.data <- marginSums(ddSX, c(2,3))
    rm(ddSX)
    gc()
    # check if any quadrature points have zero conditional intensity;
    # these do not contribute; the associated values of the sufficient
    # statistic may be Infinite and must be explicitly set to zero.
    zerocif <- (lam == 0)
    anyzerocif <- any(zerocif)
    if(anyzerocif)
      eff.data[zerocif, ] <- 0 
    ## effect of addition/deletion of U[j] on integral term in score
    changesign <- ifelse(isdata, -1, 1)
    if(!sparse) {
      ## model matrix after addition/deletion of each U[j]
      ## mombefore[i,j,] <- mom[i,]
      di <- dim(ddS)
      mombefore <- array(apply(mom, 2, rep, times=di[2]), dim=di)
      momchange <- ddS
      momchange[ , isdata, ] <- - momchange[, isdata, ]
      momafter <- mombefore + momchange
      ## effect of addition/deletion of U[j] on lambda(U[i], X)
      lamratio <- exp(tensor::tensor(momchange[,,REG,drop=FALSE], theta, 3, 1))
      lamratio <- array(lamratio, dim=dim(momafter))
      ddSintegrand <- lam * (momafter * lamratio - mombefore)
      rm(lamratio, momchange, mombefore, momafter)
      gc()
    } else {
      ## Entries are required only for pairs i,j which interact.
      ## mombefore[i,j,] <- mom[i,]
      mombefore <- mapSparseEntries(ddS, 1, mom, conform=TRUE, across=3)
      momchange <- ddS
      momchange[ , isdata, ] <- - momchange[, isdata, ]
      momafter <- mombefore + momchange
      ## lamarray[i,j,k] <- lam[i]
      lamarray <- mapSparseEntries(ddS, 1, lam, conform=TRUE, across=3)
      lamratiominus1 <- expm1(tenseur(momchange[,,REG,drop=FALSE], theta, 3, 1))
      lamratiominus1 <- as.sparse3Darray(lamratiominus1)
      ddSintegrand <- lamarray * (momafter * lamratiominus1 + momchange)
      rm(lamratiominus1, lamarray, momchange, mombefore, momafter)
      gc()
    }
    if(anyzerocif) {
      ddSintegrand[zerocif,,] <- 0
      ddSintegrand[,zerocif,] <- 0
    }
    ## integrate 
    eff.back <- changesign * tenseur(ddSintegrand, w, 1, 1)
    # total
    eff <- eff + eff.data - eff.back
    eff <- as.matrix(eff)
  } else ddSintegrand <- NULL

  # 
  if("increments" %in% what) {
    result$increm <- list(ddS=ddS,
                          ddSintegrand=ddSintegrand,
                          isdata=isdata,
                          wQ=w)
  }
  if(!any(c("leverage", "influence", "dfbetas") %in% what))
    return(result)

  # ............ compute leverage, influence, dfbetas ..............
  
  # compute basic contribution from each quadrature point
  nloc <- npoints(loc)
  switch(method,
         interpreted = {
           b <- numeric(nloc)
           for(i in seq(nloc)) {
             effi <- eff[i,, drop=FALSE]
             momi <- mom[i,, drop=FALSE]
             b[i] <- momi %*% invhess %*% t(effi)
           }
         },
         C = {
           b <- bilinearform(mom, invhess, eff)
         })
  
  # .......... leverage .............
  
  if("leverage" %in% what) {
    # values of leverage (diagonal) at points of 'loc'
    h <- b * lam
    levval <- loc %mark% h
    levval <- levval[is.finite(h)]
    levsmo <- Smooth(levval, sigma=maxnndist(loc))
    # nominal mean level
    a <- area(loc$window)
    levmean <- p/a
    lev <- list(val=levval, smo=levsmo, ave=levmean)
    result$lev <- lev
  }
  # .......... influence .............
  if("influence" %in% what) {
    # values of influence at data points
    X <- loc[isdata]
    M <- (1/p) * b[isdata]
    V <- X %mark% M
    result$infl <- V
  }
  # .......... dfbetas .............
  if("dfbetas" %in% what) {
    vex <- invhess %*% t(eff)
    switch(method,
           interpreted = {
             dis <- con <- matrix(0, nloc, ncol(mom))
             for(i in seq(nloc)) {
               vexi <- vex[,i, drop=FALSE]
               dis[i, ] <- if(isdata[i]) vexi else 0
               con[i, ] <- - lam[i] * vexi
             }
           },
           C = {
             tvex <- t(vex)
             dis <- tvex
             dis[!isdata,] <- 0
             con <- - lam  * tvex
             con[lam == 0,] <- 0
           })
    colnames(dis) <- colnames(con) <- colnames(mom)
    result$dfbetas <- msr(Q, dis[isdata, ], con)
  }
  return(result)
}

## ppmInfluenceRecombine <- function(indices, values) {
  #'  indices is a list of subset index vectors for each block
  #'  values is a list containing the computed values for each block
#           nblocks <- length(indices)
  # deriv: list(mom, score, fush, vc, hess, invhess)
  # increm: list(ddS, ddSintegrand, isdata, wQ)
  # lev: list(val, smo, ave)
  # infl: V
  # dfbetas: msr(Q, dis[isdata, ], con)
# }

## extract derivatives from covariate functions
## WARNING: these are not the score components in general

ppmDerivatives <- function(fit, what=c("gradient", "hessian"),
                            Dcovfun=NULL, loc, covfunargs=list()) {
  what <- match.arg(what)
  if(!is.null(Dcovfun)) {
    ## use provided function Dcov to compute derivatives
    Dvalues <- mpl.get.covariates(Dcovfun, loc, covfunargs=covfunargs)
    result <- as.matrix(as.data.frame(Dvalues))
    return(result)
  }
  ## any irregular parameters?
  if(length(covfunargs) == 0)
    return(NULL)
  ## Try to extract derivatives from covariate functions
  ## This often works if the functions were created by symbolic differentiation
  fvalues <- mpl.get.covariates(fit$covariates, loc, covfunargs=covfunargs,
                                need.deriv=TRUE)
  Dlist <- attr(fvalues, "derivatives")[[what]]
  if(length(Dlist) == 0)
    return(NULL)
  switch(what,
         gradient = {
           result <- do.call(cbind, unname(lapply(Dlist, as.data.frame)))
           result <- as.matrix(result)
         },
         hessian = {
           ## construct array containing Hessian matrices
           biga <- do.call(blockdiagarray, Dlist)
           ## flatten matrices 
           result <- matrix(biga, nrow=dim(biga)[1])
         })
  return(result)
}

plot.leverage.ppm <- function(x, ..., showcut=TRUE, col.cut=par("fg")) {
  fitname <- x$fitname
  defaultmain <- paste("Leverage for", fitname)
  y <- x$lev
  do.call(plot.im,
          resolve.defaults(list(y$smo),
                           list(...),
                           list(main=defaultmain)))
  if(showcut && diff(range(y$smo)) != 0) 
    do.call.matched(contour.im,
                    resolve.defaults(list(x=y$smo, levels=y$ave,
                                          add=TRUE, col=col.cut),
                                     list(...),
                                     list(drawlabels=FALSE)),
                    extrargs=c("levels", "drawlabels",
                      "labcex", "col", "lty", "lwd", "frameplot"))
  invisible(NULL)
}

plot.influence.ppm <- function(x, ...) {
  fitname <- x$fitname
  defaultmain <- paste("Influence for", fitname)
  do.call(plot.ppp,
          resolve.defaults(list(x$infl),
                           list(...),
                           list(main=defaultmain)))
}

persp.leverage.ppm <- function(x, ..., main) {
  if(missing(main)) main <- deparse(substitute(x))
  persp(as.im(x), main=main, ...)
}
  
as.im.leverage.ppm <- function(X, ...) {
  return(X$lev$smo)
}

as.ppp.influence.ppm <- function(X, ...) {
  return(X$infl)
}

as.owin.leverage.ppm <- function(W, ..., fatal=TRUE) {
  as.owin(as.im(W), ..., fatal=fatal)
}

as.owin.influence.ppm <- function(W, ..., fatal=TRUE) {
  as.owin(as.ppp(W), ..., fatal=fatal)
}

domain.leverage.ppm <- domain.influence.ppm <-
  Window.leverage.ppm <- Window.influence.ppm <-
  function(X, ...) { as.owin(X) } 

print.leverage.ppm <- function(x, ...) {
  cat("Point process leverage function\n")
  fitname <- x$fitname
  cat(paste("for model:", fitname, "\n"))
  lev <- x$lev
  cat("\nExact values:\n")
  print(lev$val)
  cat("\nSmoothed values:\n")
  print(lev$smo)
  ## for compatibility we retain the x$fit usage
  if(x$fit.is.poisson %orifnull% is.poisson(x$fit))
    cat(paste("\nAverage value:", lev$ave, "\n"))
  return(invisible(NULL))
}

print.influence.ppm <- function(x, ...) {
  cat("Point process influence measure\n")  
  fitname <- x$fitname
  cat(paste("for model:", fitname, "\n"))
  cat("\nExact values:\n")
  print(x$infl)
  return(invisible(NULL))
}

"[.leverage.ppm" <- function(x, i, ...) {
  if(missing(i)) return(x)
  y <- x$lev
  y$smo <- vi <- y$smo[i, ...]
  if(!is.im(vi)) return(vi)
  y$val <- y$val[i, ...]
  x$lev <- y
  return(x)
}

"[.influence.ppm" <- function(x, i, ...) {
  if(missing(i)) return(x)
  y <- x$infl[i, ...]
  if(!is.ppp(y)) return(y)
  x$infl <- y
  return(x)
}

shift.leverage.ppm <- function(X, ...) {
  vec <- getlastshift(shift(as.owin(X), ...))
  X$lev$val <- shift(X$lev$val, vec=vec)
  X$lev$smo <- shift(X$lev$smo, vec=vec)
  return(putlastshift(X, vec))
}

shift.influence.ppm <- function(X, ...) {
  X$infl <- shift(X$infl, ...)
  return(putlastshift(X, getlastshift(X$infl)))
}

