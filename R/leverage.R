#
#  leverage.R
#
#  leverage and influence
#
#          CORRECTED VERSION !!
#

leverage <- function(model, ...) {
  UseMethod("leverage")
}

leverage.ppm <- function(model, ...,
                         drop=FALSE, iScore=NULL, iHessian=NULL, iArgs=NULL)
{
  fitname <- short.deparse(substitute(model))
  a <- ppmInfluence(model, what="leverage", drop=drop,
                    iScore=iScore, iHessian=iHessian, iArgs=iArgs,
                    ...,
                    fitname=fitname)
  return(a$leverage)
}

influence.ppm <- function(model, ...,
                          drop=FALSE, iScore=NULL, iHessian=NULL, iArgs=NULL)
{
  fitname <- short.deparse(substitute(model))
  a <- ppmInfluence(model, what="influence", drop=drop,
                    iScore=iScore, iHessian=iHessian, iArgs=iArgs,
                    ...,
                    fitname=fitname)
  return(a$influence)
}

dfbetas.ppm <- function(model, ...,
                        drop=FALSE, iScore=NULL, iHessian=NULL, iArgs=NULL) {
  fitname <- short.deparse(substitute(model))
  a <- ppmInfluence(model, what="dfbetas", drop=drop,
                         iScore=iScore, iHessian=iHessian, iArgs=iArgs,
                     ...,
                    fitname=fitname)
  return(a$dfbetas)
}

ppmInfluence <- function(fit,
                         what=c("leverage", "influence", "dfbetas"),
                         ..., 
                         iScore=NULL, iHessian=NULL, iArgs=NULL,
                         drop=FALSE,
                         fitname=NULL) {
  stuff <- ppmInfluenceEngine(fit, what=what,
                          ..., 
                          iScore=iScore, iHessian=iHessian, iArgs=iArgs,
                          drop=drop, fitname=fitname)
  fnam <- c("fitname", "fit.is.poisson")
  result <- list()
  if("lev" %in% names(stuff)) {
    lev <- stuff[c(fnam, "lev")]
    class(lev) <- "leverage.ppm"
    result$leverage <- lev
  }
  if("infl" %in% names(stuff)) {
    infl <- stuff[c(fnam, "infl")]
    class(infl) <- "influence.ppm"
    result$influence <- infl
  }
  if(!is.null(dfb <- stuff$dfbetas)) {
    attr(dfb, "info") <- stuff[fnam]
    result$dfbetas <- dfb
  }
  other <- setdiff(names(stuff), c("lev", "infl", "dfbetas"))
  result[other] <- stuff[other]
  return(result)
}

ppmInfluenceEngine <- function(fit,
                         what=c("leverage", "influence", "dfbetas",
                           "derivatives", "increments"),
                         ...,
                         iScore=NULL, iHessian=NULL, iArgs=NULL,
                         drop=FALSE,
                         method=c("C", "interpreted"),
                         precomputed=list(),
                         sparseOK=TRUE,
                         fitname=NULL,
                         multitypeOK=FALSE,
                         entrywise = TRUE) {
  logi <- fit$method == "logi"
  if(is.null(fitname)) 
    fitname <- short.deparse(substitute(fit))
  stopifnot(is.ppm(fit))
  what <- match.arg(what, several.ok=TRUE)
  method <- match.arg(method)
  info <- list(fitname=fitname, fit.is.poisson=is.poisson(fit))
  if(is.null(iArgs))
    iArgs <- fit$covfunargs
  gotScore <- !is.null(iScore)
  gotHess <- !is.null(iHessian)
  influencecalc <- any(what %in% c("leverage", "influence", "dfbetas"))
  needHess <- gotScore && influencecalc
  if(!gotHess && needHess)
    stop("Must supply iHessian", call.=FALSE)
  sparse <- sparseOK 
  #
  # extract precomputed values if given
  theta  <- precomputed$coef   %orifnull% coef(fit)
  lam    <- precomputed$lambda %orifnull% fitted(fit, check=FALSE)
  mom    <- precomputed$mom    %orifnull% model.matrix(fit)
  # 
  p <- length(theta)
  if(logi){
    ## Intensity of dummy points
    rho <- fit$Q$param$rho %orifnull% intensity(as.ppp(fit$Q))
    logiprob <- lam / (lam + rho)
    vclist <- vcov(fit, what = "internals")
    hess <- vclist$Slog
    fush <- vclist$fisher
    invhess <- solve(hess)
    vc <- invhess %*% (vclist$Sigma1log+vclist$Sigma2log) %*% invhess
  } else{
  vc <- vcov(fit, hessian=TRUE)
  fush <- hess <- solve(vc)
  }
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

  ## evaluate additional (`irregular') components of score
  iscoremat <- ppmDerivatives(fit, "gradient", iScore, loc, covfunargs=iArgs)
  gotScore <- !is.null(iscoremat)
  needHess <- gotScore && influencecalc
  if(gotScore) {
    if(logi)
      stop("ppm influence measures are not yet implemented for method=logi with irregular parameters.")
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
               isfin <- is.finite(wlam) & matrowall(is.finite(ihessmat))
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
    if(!logi){
    hess <- fush
    invhess <- vc
    }
  }
  #
  ok <- NULL
  if(drop) {
    ok <- complete.cases(mom)
    if(all(ok)) {
      ok <- NULL
    } else {
      if((nbad <- sum(isdata[!ok])) > 0)
        warning(paste("NA value of canonical statistic at",
                      nbad, ngettext(nbad, "data point", "data points")),
                call.=FALSE)
      Q <- Q[ok]
      mom <- mom[ok, , drop=FALSE]
      loc <- loc[ok]
      lam <- lam[ok]
      w   <- w[ok]
      isdata <- isdata[ok]
    }
  } 

  # ........  start assembling results .....................
  # 
  result <- as.list(info)
  class(result) <- "ppmInfluence"
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
  #  ........ Poisson case ..................................
  eff <- mom
  # ........  Gibbs case ....................................
  ## second order interaction terms
  ddS <- ddSintegrand <- NULL
  if(!is.poisson(fit)) {
    ## goal is to compute these effect matrices:
    eff.data <- eff.back  <- matrix(0, nrow(eff), ncol(eff),
                                    dimnames=dimnames(eff))
    U <- union.quad(Q)
    nU <- npoints(U)
    zerocif <- (lam == 0)
    anyzerocif <- any(zerocif)
    ## decide whether to split into blocks
    nX <- Q$data$n
    nD <- Q$dummy$n
    bls <- quadBlockSizes(nX, nD, p, announce=TRUE)
    nblocks    <- bls$nblocks
    nperblock  <- bls$nperblock
    ##
    if(nblocks > 1 && ("increments" %in% what)) {
      warning("Oversize quadrature scheme: cannot return array of increments",
              call.=FALSE)
      what <- setdiff(what, "increments")
    }
    R <- reach(fit)
    ## indices into original quadrature scheme
    whichok <- if(!is.null(ok)) which(ok) else seq_len(nX+nD) 
    whichokdata <- whichok[isdata]
    whichokdummy <- whichok[!isdata]
    ## loop 
    for(iblock in 1:nblocks) {
      first <- min(nD, (iblock - 1) * nperblock + 1)
      last  <- min(nD, iblock * nperblock)
      # corresponding subset of original quadrature scheme
      if(!is.null(ok) || nblocks > 1) {
        ## subset for which we will compute the effect
        current <- c(whichokdata, whichokdummy[first:last])
        ## find neighbours, needed for calculation
        other <- setdiff(seq_len(nU), current)
        crx <- crosspairs(U[current], U[other], R, what="indices")
        nabers <- other[unique(crx$j)]
        ## subset actually requested
        requested <- c(current, nabers)
        ## corresponding stuff ('B' for block)
        isdataB <- isdata[requested]
        zerocifB <- zerocif[requested]
        anyzerocifB <- any(zerocifB)
        momB <- mom[requested, , drop=FALSE]
        lamB <- lam[requested]
        if(logi) logiprobB <- logiprob[requested]
        wB <- w[requested]
        currentB <- seq_along(current)
      } else {
        requested <- NULL
        isdataB <- isdata
        zerocifB <- zerocif
        anyzerocifB <- anyzerocif
        momB <- mom
        lamB <- lam
        if(logi) logiprobB <- logiprob
        wB <- w
      }
      ## compute second order terms 
      ## ddS[i,j, ] = Delta_i Delta_j S(x)
      ddS <- deltasuffstat(fit, restrict = "first", dataonly=FALSE,
                           quadsub=requested, sparseOK=sparse)
      ## 
      if(is.null(ddS)) {
        warning("Second order interaction terms are not implemented",
                " for this model; they are treated as zero", call.=FALSE)
        break
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
      ## effect of addition/deletion of U[j]
      ## on score contribution from data points (sum automatically restricted to
      ## interior for border correction by earlier call to
      ## deltasuffstat(..., restrict = "first"))
      ddSX <- ddS[isdataB, , , drop=FALSE]
      eff.data.B <- marginSums(ddSX, c(2,3))
      ## check if any quadrature points have zero conditional intensity;
      ## these do not contribute; the associated values of the sufficient
      ## statistic may be Infinite and must be explicitly set to zero.
      if(anyzerocifB)
        eff.data.B[zerocifB, ] <- 0
      ## save results for current subset of quadrature points 
      if(is.null(requested)) {
        eff.data <- eff.data.B
      } else {
        eff.data[current,] <- as.matrix(eff.data.B[currentB,,drop=FALSE])
      }
      ## 
      rm(ddSX, eff.data.B)
      ## effect of addition/deletion of U[j] on integral term in score
      changesignB <- ifelse(isdataB, -1, 1)
      if(!sparse) {
        if(logi){
          stop("Non-sparse method is not implemented for method = 'logi'!")
        } else{
          ## model matrix after addition/deletion of each U[j]
          ## mombefore[i,j,] <- mom[i,]
          di <- dim(ddS)
          mombefore <- array(apply(momB, 2, rep, times=di[2]), dim=di)
          momchange <- ddS
          momchange[ , isdataB, ] <- - momchange[, isdataB, ]
          momafter <- mombefore + momchange
          ## effect of addition/deletion of U[j] on lambda(U[i], X)
          if(gotScore){
            lamratio <- exp(tensor::tensor(momchange[,,REG,drop=FALSE],
                                           theta, 3, 1))
          } else{
            lamratio <- exp(tensor::tensor(momchange, theta, 3, 1))
          }
          lamratio <- array(lamratio, dim=dim(momafter))
          ddSintegrand <- lamB * (momafter * lamratio - mombefore)
          rm(lamratio)
        }
        rm(momchange, mombefore, momafter)
        gc()
      } else {
        if(logi){
          ## Delta suff. stat. with sign change for data/dummy (sparse3Darray)
          momchange <- ddS
          momchange[ , isdataB, ] <- - momchange[, isdataB, ]
          ## theta^T %*% ddS (with sign -1/+1 according to data/dummy) as triplet sparse matrix
          if(gotScore){
            momchangeeffect <- tenseur(momchange[,,REG,drop=FALSE], theta, 3, 1)
          } else{
            momchangeeffect <- tenseur(momchange, theta, 3, 1)
          }
          momchangeeffect <- expandSparse(momchangeeffect, n = dim(ddS)[3], across = 3)
          ijk <- SparseIndices(momchangeeffect)
          ## Entrywise calculations below
          momchange <- as.numeric(momchange[ijk])
          ## Transform to change in probability
          expchange <- exp(momchangeeffect$x)
          lamBi <- lamB[ijk$i]
          logiprobBi <- logiprobB[ijk$i]
          changesignBj <- changesignB[ijk$j]
          pchange <- changesignBj*(lamBi * expchange / (lamBi*expchange + rho) - logiprobBi)
          mombefore <- mom[cbind(ijk$i,ijk$k)]
          ## Note: changesignBj * momchange == as.numeric(ddS[ijk])
          ddSintegrand <- (mombefore + momchange) * pchange + logiprobBi * changesignBj * momchange
          ddSintegrand <- sparse3Darray(i = ijk$i, j = ijk$j, k = ijk$k, x = ddSintegrand,
                                        dims = dim(ddS))
        } else{
          if(entrywise){
            momchange <- ddS
            momchange[ , isdataB, ] <- - momchange[, isdataB, ]
            if(gotScore){
              lamratiominus1 <- expm1(tenseur(momchange[,,REG,drop=FALSE],
                                              theta, 3, 1))
            } else{
              lamratiominus1 <- expm1(tenseur(momchange, theta, 3, 1))
            }
            lamratiominus1 <- expandSparse(lamratiominus1, n = dim(ddS)[3], across = 3)
            ijk <- SparseIndices(lamratiominus1)
            ## Everything entrywise with ijk now:
            # lamratiominus1 <- lamratiominus1[cbind(ijk$i, ijk$j)]
            lamratiominus1 <- as.numeric(lamratiominus1$x)
            momchange <- as.numeric(momchange[ijk])
            mombefore <- momB[cbind(ijk$i, ijk$k)]
            momafter <- mombefore + momchange
            ## lamarray[i,j,k] <- lam[i]
            lamarray <- lamB[ijk$i]
            ddSintegrand <- lamarray * (momafter * lamratiominus1 + momchange)
            ddSintegrand <- sparse3Darray(i = ijk$i, j = ijk$j, k = ijk$k, x = ddSintegrand,
                                          dims = dim(ddS))
          } else{
            ## Entries are required only for pairs i,j which interact.
            ## mombefore[i,j,] <- mom[i,]
            mombefore <- mapSparseEntries(ddS, 1, momB, conform=TRUE, across=3)
            momchange <- ddS
            momchange[ , isdataB, ] <- - momchange[, isdataB, ]
            momafter <- evalSparse3Dentrywise(mombefore + momchange)
            ## lamarray[i,j,k] <- lam[i]
            lamarray <- mapSparseEntries(ddS, 1, lamB, conform=TRUE, across=3)
            if(gotScore){
              lamratiominus1 <- expm1(tenseur(momchange[,,REG,drop=FALSE],
                                              theta, 3, 1))
            } else{
              lamratiominus1 <- expm1(tenseur(momchange,theta, 3, 1))
            }
            lamratiominus1 <- expandSparse(lamratiominus1, n = dim(ddS)[3], across = 3)
            ddSintegrand <- evalSparse3Dentrywise(lamarray * (momafter* lamratiominus1 + momchange))
          }
          rm(lamratiominus1, lamarray, momafter)
        }
        rm(momchange, mombefore)
      }
      if(anyzerocifB) {
        ddSintegrand[zerocifB,,] <- 0
        ddSintegrand[,zerocifB,] <- 0
      }
      ## integrate
      if(logi){
        # eff.back.B <- tenseur(ddSintegrand, rep(1, length(wB)), 1, 1)
        eff.back.B <- marginSums(ddSintegrand, c(2,3))
      } else{
        eff.back.B <- changesignB * tenseur(ddSintegrand, wB, 1, 1)
      }
      ## save contribution
      if(is.null(requested)) {
        eff.back <- eff.back.B
      } else {
        eff.back[current,] <- as.matrix(eff.back.B[currentB, , drop=FALSE])
      }
    }
    
    ## total
    eff <- eff + eff.data - eff.back
    eff <- as.matrix(eff)
  }
  
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
    ## values of leverage (diagonal) at points of 'loc'
    h <- b * lam
    ok <- is.finite(h)
    if(mt <- is.multitype(loc))
      h <- data.frame(leverage=h, type=marks(loc))
    levval <- (loc %mark% h)[ok]
    levvaldum <- levval[!isdata[ok]]
    if(!mt) {
      levsmo <- Smooth(levvaldum, sigma=maxnndist(loc))
    } else {
      levsplitdum <- split(levvaldum, reduce=TRUE)
      levsmo <- Smooth(levsplitdum, sigma=max(sapply(levsplitdum, maxnndist)))
    }
    ## nominal mean level
    a <- area(Window(loc)) * markspace.integral(loc)
    levmean <- p/a
    lev <- list(val=levval, smo=levsmo, ave=levmean)
    result$lev <- lev
  }
  # .......... influence .............
  if("influence" %in% what) {
    if(logi){
      X <- loc
      effX <- as.numeric(isdata) * eff - mom * logiprob
    } else{
      # values of influence at data points
      X <- loc[isdata]
      effX <- eff[isdata, ,drop=FALSE]
    }
    M <- (1/p) * quadform(effX, invhess)
    if(logi || is.multitype(X)) {
      # result will have several columns of marks
      M <- data.frame(influence=M)
      if(logi) M$isdata <- factor(isdata, levels = c(TRUE, FALSE), labels = c("data", "dummy"))
      if(is.multitype(X)) M$type <- marks(X)
    } 
    V <- X %mark% M
    result$infl <- V
  }
  # .......... dfbetas .............
  if("dfbetas" %in% what) {
    if(logi){
      M <- as.numeric(isdata) * eff - mom * logiprob
      M <- t(invhess %*% t(M))
      Mdum <- M
      Mdum[isdata,] <- 0
      Mdum <- Mdum / w.quad(Q)
      result$dfbetas <- msr(Q, M[isdata, ], Mdum)
    } else{
      vex <- invhess %*% t(mom)
      dex <- invhess %*% t(eff)
      switch(method,
             interpreted = {
               dis <- con <- matrix(0, nloc, ncol(mom))
               for(i in seq(nloc)) {
                 vexi <- vex[,i, drop=FALSE]
                 dexi <- dex[,i, drop=FALSE]
                 dis[i, ] <- if(isdata[i]) dexi else 0
                 con[i, ] <- - lam[i] * vexi
               }
             },
             C = {
               dis <- t(dex)
               dis[!isdata,] <- 0
               con <- - lam * t(vex)
               con[lam == 0,] <- 0
             })
      colnames(dis) <- colnames(con) <- colnames(mom)
      # result is a vector valued measure
      result$dfbetas <- msr(Q, dis[isdata, ], con)
    }
  }
  return(result)
}

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
           result <- matrix(biga, nrow=dim(biga)[1L])
         })
  return(result)
}

plot.leverage.ppm <- local({

  plot.leverage.ppm <- function(x, ..., showcut=TRUE, col.cut=par("fg"),
                                multiplot=TRUE) {
    fitname <- x$fitname
    defaultmain <- paste("Leverage for", fitname)
    y <- x$lev
    smo <- y$smo
    ave <- y$ave
    if(!multiplot && inherits(smo, "imlist")) {
      ave <- ave * length(smo)
      smo <- Reduce("+", smo)
      defaultmain <- c(defaultmain, "(sum over all types of point)")
    }
    if(is.im(smo)) {
      do.call(plot.im,
              resolve.defaults(list(smo),
                               list(...),
                               list(main=defaultmain)))
      if(showcut)
        onecontour(0, x=smo, ..., level.cut=ave, col.cut=col.cut)
    } else if(inherits(smo, "imlist")) {
      xtra <- list(panel.end=onecontour,
                   panel.end.args=list(level.cut=ave, col.cut=col.cut))
      do.call(plot.solist,
              resolve.defaults(list(smo),
                               list(...),
                               list(main=defaultmain),
                               if(showcut) xtra else list()))
    } 
    invisible(NULL)
  }

  onecontour <- function(i, x, ..., level.cut, col.cut) {
    if(diff(range(x)) > 0)
      do.call.matched(contour.im,
                      resolve.defaults(list(x=x, levels=level.cut,
                                            add=TRUE, col=col.cut),
                                       list(...),
                                       list(drawlabels=FALSE)),
                      extrargs=c("levels", "drawlabels",
                        "labcex", "col", "lty", "lwd", "frameplot"))
  }

  plot.leverage.ppm
})
                           
plot.influence.ppm <- function(x, ..., multiplot=TRUE) {
  fitname <- x$fitname
  defaultmain <- paste("Influence for", fitname)
  y <- x$infl
  if(multiplot && isTRUE(ncol(marks(y)) > 1)) {
    # apart from the influence value, there may be additional columns of marks
    # containing factors: {type of point}, { data vs dummy in logistic case }
    ma <- as.data.frame(marks(y))
    fax <- sapply(ma, is.factor)
    nfax <- sum(fax)
    if(nfax == 1) {
      # split on first available factor, and remove this factor
      y <- split(y, reduce=TRUE)
    } else if(nfax > 1) {
      # several factors: split according to them all, and remove them all
      f.all <- do.call(interaction, ma[fax])
      z <- y %mark% ma[,!fax]
      y <- split(z, f.all)
    }
  }
  do.call(plot,
          resolve.defaults(list(y),
                           list(...),
                           list(main=defaultmain,
                                multiplot=multiplot,
                                which.marks=1)))
}

persp.leverage.ppm <- function(x, ..., main) {
  if(missing(main)) main <- deparse(substitute(x))
  y <- as.im(x)
  if(is.im(y)) return(persp(y, main=main, ...))
  pa <- par(ask=TRUE)
  lapply(y, persp, main=main, ...)
  par(pa)
  return(invisible(NULL))
}
  
as.im.leverage.ppm <- function(X, ...) {
  return(X$lev$smo) # could be either an image or a list of images
}

as.function.leverage.ppm <- function(x, ...) {
  X <- x$lev$val
  S <- ssf(unmark(X), marks(X))
  return(as.function(S))
}

as.ppp.influence.ppm <- function(X, ...) {
  return(X$infl)
}

as.owin.leverage.ppm <- function(W, ..., fatal=TRUE) {
  y <- as.im(W)
  if(inherits(y, "imlist")) y <- y[[1L]]
  as.owin(y, ..., fatal=fatal)
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
  smoi <- if(is.im(y$smo)) y$smo[i, ...] else solapply(y$smo, "[", i=i, ...)
  if(!inherits(smoi, c("im", "imlist"))) return(smoi)
  y$smo <- smoi
  y$val <- y$val[i, ...]
  y$ave <- if(is.im(smoi)) mean(smoi) else mean(sapply(smoi, mean))
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
  smo <- X$lev$smo
  X$lev$smo <-
    if(is.im(smo)) shift(smo, vec=vec) else solapply(smo, shift, vec=vec)
  return(putlastshift(X, vec))
}

shift.influence.ppm <- function(X, ...) {
  X$infl <- shift(X$infl, ...)
  return(putlastshift(X, getlastshift(X$infl)))
}

