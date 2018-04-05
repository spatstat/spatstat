#
#  leverage.R
#
#  leverage and influence
#
#  $Revision: 1.109 $ $Date: 2018/04/05 03:30:11 $
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
  class(result) <- "ppmInfluence"
  return(result)
}

leverage.ppmInfluence <- function(model, ...) { model$leverage }
influence.ppmInfluence <- function(model, ...) { model$influence }
dfbetas.ppmInfluence <- function(model, ...) { model$dfbetas }

avenndist <- function(X) mean(nndist(X))

## ...............  main workhorse ....................................

ppmInfluenceEngine <- function(fit,
                         what=c("leverage", "influence", "dfbetas",
                           "score", "derivatives", "increments"),
                         ...,
                         iScore=NULL, iHessian=NULL, iArgs=NULL,
                         drop=FALSE,
                         method=c("C", "interpreted"),
                         fine=FALSE,
                         precomputed=list(),
                         sparseOK=TRUE,
                         fitname=NULL,
                         multitypeOK=FALSE,
                         entrywise = TRUE,
                         matrix.action = c("warn", "fatal", "silent"),
                         dimyx=NULL, eps=NULL,
                         geomsmooth = TRUE) {
  if(is.null(fitname)) 
    fitname <- short.deparse(substitute(fit))

  ## type of calculation to be performed
  method <- match.arg(method)
  what <- match.arg(what, several.ok=TRUE)
  matrix.action <- match.arg(matrix.action)

  influencecalc <- any(what %in% c("leverage", "influence", "dfbetas"))
  hesscalc <- influencecalc || any(what == "derivatives")
  sparse <- sparseOK 
  target <- paste(what, collapse=",")
  
  ## ...........  collect information about the model .................
  
  stopifnot(is.ppm(fit))

  #' ensure object contains GLM fit
  if(!hasglmfit(fit)) {
    fit <- update(fit, forcefit=TRUE)
    precomputed <- list()
  }

  #' type of interpoint interaction
  fit.is.poisson <- is.poisson(fit)
  hasInf <- !fit.is.poisson && !identical(fit$interaction$hasInf, FALSE)

  #' estimating function 
  fitmethod <- fit$method
  logi <- (fitmethod == "logi")
  pseudo <- (fitmethod == "mpl")
  if(!logi && !pseudo) {
    warning(paste("Model was fitted with method =", dQuote(fitmethod),
                  "but is treated as having been fitted by maximum",
		  if(fit.is.poisson) "likelihood" else "pseudolikelihood",
		  "for leverage/influence calculation"),
	    call.=FALSE)
    pseudo <- TRUE
  }

  ## Detect presence of irregular parameters
  if(is.null(iArgs))
    iArgs <- fit$covfunargs
  gotScore <- !is.null(iScore)
  gotHess <- !is.null(iHessian)
  needHess <- gotScore && hesscalc  # may be updated later
  if(!gotHess && needHess)
    stop("Must supply iHessian", call.=FALSE)

  #' ...................  evaluate basic terms ....................
  
  ## extract values from model, using precomputed values if given
  theta  <- precomputed$coef   %orifnull% coef(fit)
  lampos <- precomputed$lambda %orifnull% fitted(fit, ignore.hardcore=hasInf,
                                                      check=FALSE)
  mom    <- precomputed$mom    %orifnull% model.matrix(fit, splitInf=hasInf)

  ## 'lampos' is positive part of cif
  ## 'lam' is full model cif including zeroes
  lam <- lampos
  zerocif <- attr(mom, "-Inf") %orifnull% logical(nrow(mom))
  anyzerocif <- any(zerocif)
  if(hasInf && anyzerocif)
    lam[zerocif] <- 0
  
  p <- length(theta)
  Q <- quad.ppm(fit)
  w <- w.quad(Q)
  loc <- union.quad(Q)
  isdata <- is.data(Q)
  mt <- is.multitype(loc)
  if(length(w) != length(lam))
    stop(paste("Internal error: length(w) = ", length(w),
               "!=", length(lam), "= length(lam)"),
         call.=FALSE)

  ## smoothing bandwidth and resolution for smoothed images of densities
  smallsigma <- if(!mt) avenndist(loc) else max(sapply(split(loc), avenndist))
  ## previously used 'maxnndist' instead of 'avenndist'
  if(is.null(dimyx) && is.null(eps)) 
    eps <- sqrt(prod(sidelengths(Frame(loc))))/256

  #' ...............  evaluate Hessian of regular parameters ................

  ## domain of composite likelihood
  ## (e.g. eroded window in border correction)
  inside <- getglmsubset(fit) %orifnull% rep(TRUE, npoints(loc))
  
  ## extract negative Hessian matrix of regular part of log composite likelihood
  ##  hess = negative Hessian H
  ##  fgrad = Fisher-scoring-like gradient G = estimate of E[H]

  if(logi) {
    ## ..............    logistic composite likelihood ......................
    ## Intensity of dummy points
    rho <- fit$Q$param$rho %orifnull% intensity(as.ppp(fit$Q))
    logiprob <- lampos / (lampos + rho)
    vclist <- vcov(fit, what = "internals", fine=fine, matrix.action="silent")
    hess <- vclist$Slog
    fgrad <- vclist$fisher
    invhess <- if(is.null(hess)) NULL else checksolve(hess, "silent")
    invfgrad <- if(is.null(fgrad)) NULL else checksolve(fgrad, "silent")
    if(is.null(invhess) || is.null(invfgrad)) {
      #' use more expensive estimate of variance terms
      vclist <- vcov(fit, what = "internals", fine=TRUE,
                     matrix.action=matrix.action)
      hess <- vclist$Slog
      fgrad <- vclist$fisher
      #' try again - exit if really singular
      invhess <- checksolve(hess, matrix.action, "Hessian", target)
      invfgrad <- checksolve(fgrad, matrix.action, "gradient matrix", target)
    }
#    vc <- invhess %*% (vclist$Sigma1log+vclist$Sigma2log) %*% invhess
  } else {
    ## ..............  likelihood or pseudolikelihood ....................
    invfgrad <- vcov(fit, hessian=TRUE, fine=fine, matrix.action="silent")
    fgrad <- hess <- if(is.null(invfgrad) || anyNA(invfgrad)) NULL else
                     checksolve(invfgrad, "silent")
    if(is.null(fgrad)) {
      invfgrad <- vcov(fit, hessian=TRUE, fine=TRUE,
                       matrix.action=matrix.action)
      fgrad <- hess <- checksolve(invfgrad, matrix.action, "Hessian", target)
    }
  }

  #' ...............  augment Hessian  ...................

  ## evaluate additional (`irregular') components of score, if any
  iscoremat <- ppmDerivatives(fit, "gradient", iScore, loc, covfunargs=iArgs)
  gotScore <- !is.null(iscoremat)
  needHess <- gotScore && hesscalc
  if(!gotScore) {
    REG <- 1:ncol(mom)
  } else {
    ## count regular and irregular parameters
    nreg <- ncol(mom)
    nirr <- ncol(iscoremat)
    ## add extra columns to model matrix
    mom <- cbind(mom, iscoremat)
    REG <- 1:nreg
    IRR <- nreg + 1:nirr
    ## evaluate additional (`irregular') entries of Hessian
    ihessmat <- if(!needHess) NULL else
                ppmDerivatives(fit, "hessian", iHessian, loc, covfunargs=iArgs)
    if(gotHess <- !is.null(ihessmat)) {
      ## recompute negative Hessian of log PL and its mean
      fgrad <- hessextra <- matrix(0, ncol(mom), ncol(mom))
    }  
    if(pseudo) {
      ## ..............  likelihood or pseudolikelihood ....................
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
                     fgrad <- fgrad + v1
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
                 hess <- fgrad + hessextra
                 invhess <- checksolve(hess, matrix.action, "Hessian", target)
               } else {
                 invhess <- hess <- NULL
               }
             },
             C = {
               wlam <- lam * w
               fgrad <- sumouter(mom, wlam)
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
                 hess <- fgrad + hessextra
                 invhess <- checksolve(hess, matrix.action, "Hessian", target)
               } else {
                 invhess <- hess <- NULL
               }
             })
    } else {
      ## ..............  logistic composite likelihood ....................
      switch(method,
             interpreted = {
	       oweight <- logiprob * (1 - logiprob)
	       hweight <- ifelse(isdata, -(1 - logiprob), logiprob)
               for(i in seq(loc$n)) {
                 ## outer product of score 
                 momi <- mom[i, ]
                 v1 <- outer(momi, momi, "*") * oweight[i]
                 if(all(is.finite(v1)))
                   fgrad <- fgrad + v1
		 ## Hessian term
                 ## contributions nonzero for irregular parameters
                 if(gotHess) {
                   v2 <- hweight[i] *
		         matrix(as.numeric(ihessmat[i,]), nirr, nirr)
                   if(all(is.finite(v2)))
                     hessextra[IRR, IRR] <- hessextra[IRR, IRR] + v2
                 }
               }
	       if(gotHess) {
                 hess <- fgrad + hessextra
                 invhess <- checksolve(hess, matrix.action, "Hessian", target)
               } else {
                 invhess <- hess <- NULL
	       }
             },
             C = {
	       oweight <- logiprob * (1 - logiprob)
	       hweight <- ifelse(isdata, -(1 - logiprob), logiprob)
               fgrad <- sumouter(mom, oweight)
               if(gotHess) {
                 # Hessian term
                 isfin <- is.finite(hweight) & matrowall(is.finite(ihessmat))
                 vcontrib <-
                   if(all(isfin)) hweight %*% ihessmat else
                               hweight[isfin] %*% ihessmat[isfin,, drop=FALSE]
                 hessextra[IRR, IRR] <-
                   hessextra[IRR, IRR] + matrix(vcontrib, nirr, nirr)
                 hess <- fgrad + hessextra
                 invhess <- checksolve(hess, matrix.action, "Hessian", target)
               } else {
                 invhess <- hess <- NULL
               }
             })
    }
    invfgrad <- checksolve(fgrad, matrix.action, "gradient matrix", target)
  }
  
  if(!needHess) {
    if(pseudo){
    hess <- fgrad
    invhess <- invfgrad
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
      inside <- inside[ok]
    }
  } 

  # ........  start assembling result .....................
  result <- list(fitname=fitname, fit.is.poisson=fit.is.poisson)

  if(any(c("score", "derivatives") %in% what)) {
    ## calculate the composite score
    rawmean <- if(logi) logiprob else (lam * w)
    rawresid <- isdata - rawmean
    score <- matrix(rawresid, nrow=1) %*% mom

    if("score" %in% what)
      result$score <- score
    if("derivatives" %in% what) 
      result$deriv <- list(mom=mom, score=score,
                           fgrad=fgrad, invfgrad=invfgrad,
                           hess=hess, invhess=invhess)
    if(all(what %in% c("score", "derivatives")))
      return(result)
  }


  ## :::::::::::::::  compute second order terms :::::::::::::


  ##  >>>  set model matrix to zero outside the domain <<<
  mom[!inside, ] <- 0
  
  ## compute effect of adding/deleting each quadrature point
  if(fit.is.poisson) {
    ##  ........ Poisson case ..................................
    eff <- mom
    ddS <- ddSintegrand <- NULL
  } else {
    ## ........  Gibbs case ....................................
    ## initialise
    eff <- mom
    ## second order interaction terms
    ##    columns index the point being added/deleted
    ##    rows index the points affected
    ## goal is to compute these effect matrices:
    eff.data <- eff.back  <- matrix(0, nrow(eff), ncol(eff),
                                    dimnames=dimnames(eff))
    ## 
    U <- union.quad(Q)
    nU <- npoints(U)
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
    
    ## {{{{{{{{{{{{{   L O O P   }}}}}}}}}}}}}
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
        changesignB <- ifelse(isdataB, -1, 1)
        zerocifB <- zerocif[requested]
        anyzerocifB <- any(zerocifB)
        momB <- mom[requested, , drop=FALSE]
        lamB <- lam[requested]
        #' unused:
        #'     insideB <- inside[requested]
        #'     lamposB <- lampos[requested]
        if(logi) logiprobB <- logiprob[requested]
        wB <- w[requested]
        currentB <- seq_along(current)
      } else {
        requested <- NULL
        isdataB <- isdata
        changesignB <- ifelse(isdataB, -1, 1)
        zerocifB <- zerocif
        anyzerocifB <- anyzerocif
        momB <- mom
        lamB <- lam
        #'  unused:
        #'     insideB <- inside
	#'     lamposB <- lampos
        if(logi) logiprobB <- logiprob
        wB <- w
      }
      ## compute second order terms 
      ## ddS[i,j, ] = Delta_i Delta_j S(x)
      ddS <- deltasuffstat(fit, restrict = "first", dataonly=FALSE,
                           quadsub=requested, sparseOK=sparse,
			   splitInf=hasInf,
                           force=TRUE, warn.forced=TRUE)
      ## 
      if(is.null(ddS)) {
        warning("Second order interaction terms are not implemented",
                " for this model; they are treated as zero", call.=FALSE)
        break
      } else {
        sparse <- inherits(ddS, "sparse3Darray")
	if(hasInf) {
          deltaInf <- attr(ddS, "deltaInf")
          hasInf <- !is.null(deltaInf)
          if(hasInf) sparse <- sparse && inherits(deltaInf, "sparseMatrix")
	}
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
      ## ^^^^^^^^^^^^^^^^^  second term in DeltaScore ^^^^^^^^^^^^^^^^^^^^
      ## effect of addition/deletion of U[j]
      ## on score contribution from data points (sum automatically restricted to
      ## interior for border correction by earlier call to
      ## deltasuffstat(..., restrict = "first"))
      ddSX <- ddS[isdataB, , , drop=FALSE]
      eff.data.B <- marginSums(ddSX, c(2,3))
      ## check if any quadrature points have zero conditional intensity;
      ## these do not contribute to this term
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
      ## ^^^^^^^^^^^^^^^^^  third term in DeltaScore ^^^^^^^^^^^^^^^^^^^^
      ## effect of addition/deletion of U[j] on integral term in score
      if(!sparse) {
        ## ::::::::::::::: full arrays, simpler code :::::::::::::::::::
        if(pseudo) {
          ## ---------------  likelihood or pseudolikelihood -----------
          ## model matrix after addition/deletion of each U[j]
          ## mombefore[i,j,] <- mom[i,]
          di <- dim(ddS)
          mombefore <- array(apply(momB, 2, rep, times=di[2]), dim=di)
          momchange <- ddS
          momchange[ , isdataB, ] <- - momchange[, isdataB, ]
          momafter <- mombefore + momchange
          ## effect of addition/deletion of U[j] on lambda(U[i], X)
          if(gotScore) {
            lamratio <- exp(tensor::tensor(momchange[,,REG,drop=FALSE],
                                           theta, 3, 1))
          } else {
            lamratio <- exp(tensor::tensor(momchange, theta, 3, 1))
          }
          lamratio <- array(lamratio, dim=dim(momafter))
          if(!hasInf) {
            #' cif is positive
            ddSintegrand <- lamB * (momafter * lamratio - mombefore)
          } else {
            #' cif can be zero
            zerobefore <- matrix(zerocifB, di[1], di[2])
            zerochange <- deltaInf * 1L
            zerochange[ , isdataB] <- - zerochange[ , isdataB]
            zeroafter  <- zerobefore + zerochange
            momZbefore <- mombefore
            momZbefore[ , zerocifB, ] <- 0
            IJK <- unname(which(array(zeroafter == 1, dim=di), arr.ind=TRUE))
            momZafter <- momafter
            momZafter[IJK] <- 0
            momZchange <- momZafter- momZbefore
            ddSintegrand <- lamB * (momZafter * lamratio - momZbefore)
          }
          rm(momchange, mombefore, momafter, lamratio)
        } else {
          ## ---------------  logistic composite likelihood -----------
          stop("Non-sparse method is not implemented for method = 'logi'!")
        }
        gc()
      } else {
        ## ::::::::::::::::::   sparse arrays   ::::::::::::::::::::::::
        if(logi) {
          ## ---------------  logistic composite likelihood -----------
          ## Delta suff. stat. with sign change for data/dummy (sparse3Darray)
          momchange <- ddS
          momchange[ , isdataB, ] <- - momchange[, isdataB, ]
          ## Evaluate theta^T %*% ddS (with sign -1/+1 according to data/dummy)
          ## as triplet sparse matrix
          if(gotScore){
            momchangeeffect <- tenseur(momchange[,,REG,drop=FALSE], theta, 3, 1)
          } else{
            momchangeeffect <- tenseur(momchange, theta, 3, 1)
          }
          ## Copy to each slice 
          momchangeeffect <- expandSparse(momchangeeffect,
                                          n = dim(ddS)[3], across = 3)
          ijk <- SparseIndices(momchangeeffect)
          ## Entrywise calculations below
          momchange <- as.numeric(momchange[ijk])
          mombefore <- mom[cbind(ijk$i,ijk$k)]
          momafter  <- mombefore + momchange
          ## Transform to change in probability
          expchange <- exp(momchangeeffect$x)
          lamBi <- lamB[ijk$i]
          logiprobBi <- logiprobB[ijk$i]
          changesignBj <- changesignB[ijk$j]
          pchange <- changesignBj*(lamBi * expchange / (lamBi*expchange + rho) - logiprobBi)
          ## Note: changesignBj * momchange == as.numeric(ddS[ijk])
          if(!hasInf) {
            #' cif is positive
            ddSintegrand <- momafter * pchange +
              logiprobBi * changesignBj * momchange
          } else {
            #' cif can be zero
            isdataBj <- isdataB[ijk$j]
            zerobefore <- as.logical(zerocifB[ijk$i])
            zerochange <- as.logical(deltaInf[cbind(ijk$i, ijk$j)])
            zerochange[isdataBj] <- - zerochange[isdataBj]
            zeroafter <- zerobefore + zerochange
            momZbefore <- ifelse(zerobefore, 0, mombefore)
            momZafter  <- ifelse(zeroafter,  0, momafter)
            momZchange <- momZafter - momZbefore
            ddSintegrand <- momZafter * pchange +
              logiprobBi * changesignBj * momZchange
          }
          ddSintegrand <- sparse3Darray(i = ijk$i, j = ijk$j, k = ijk$k,
                                        x = ddSintegrand,
                                        dims = dim(ddS))
        } else{
          ## ---------------  likelihood or pseudolikelihood -----------
          if(entrywise) {
            ## ......  sparse arrays, using explicit indices ......
            momchange <- ddS
            momchange[ , isdataB, ] <- - momchange[, isdataB, ]
            if(gotScore){
              lamratiominus1 <- expm1(tenseur(momchange[,,REG,drop=FALSE],
                                              theta, 3, 1))
            } else{
              lamratiominus1 <- expm1(tenseur(momchange, theta, 3, 1))
            }
            lamratiominus1 <- expandSparse(lamratiominus1,
                                           n = dim(ddS)[3], across = 3)
            ijk <- SparseIndices(lamratiominus1)
            ## Everything entrywise with ijk now:
            # lamratiominus1 <- lamratiominus1[cbind(ijk$i, ijk$j)]
            lamratiominus1 <- as.numeric(lamratiominus1$x)
            momchange <- as.numeric(momchange[ijk])
            mombefore <- momB[cbind(ijk$i, ijk$k)]
            momafter <- mombefore + momchange
            ## lamarray[i,j,k] <- lam[i]
            lamarray <- lamB[ijk$i]
            if(!hasInf) {
              #' cif is positive
              ddSintegrand <- lamarray * (momafter * lamratiominus1 + momchange)
            } else {
              #' cif can be zero
              isdataBj <- isdataB[ijk$j]
              zerobefore <- as.logical(zerocifB[ijk$i])
              zerochange <- as.logical(deltaInf[cbind(ijk$i, ijk$j)])
              zerochange[isdataBj] <- - zerochange[isdataBj]
              zeroafter <- zerobefore + zerochange
              momZbefore <- ifelse(zerobefore, 0, mombefore)
              momZafter  <- ifelse(zeroafter,  0, momafter)
              momZchange <- momZafter - momZbefore
              ddSintegrand <- lamarray*(momZafter*lamratiominus1 + momZchange)
            }
            ddSintegrand <- sparse3Darray(i = ijk$i, j = ijk$j, k = ijk$k,
                                          x = ddSintegrand,
                                          dims = dim(ddS))
          } else {
            ## ......  sparse array code ......
            ## Entries are required only for pairs i,j which interact.
            ## mombefore[i,j,] <- mom[i,]
            mombefore <- mapSparseEntries(ddS, 1, momB, conform=TRUE, across=3)
            momchange <- ddS
            momchange[ , isdataB, ] <- - momchange[, isdataB, ]
            ##            momafter <- evalSparse3Dentrywise(mombefore + momchange)
            momafter <- mombefore + momchange
            ## lamarray[i,j,k] <- lam[i]
            lamarray <- mapSparseEntries(ddS, 1, lamB, conform=TRUE, across=3)
            if(gotScore){
              lamratiominus1 <- expm1(tenseur(momchange[,,REG,drop=FALSE],
                                              theta, 3, 1))
            } else{
              lamratiominus1 <- expm1(tenseur(momchange,theta, 3, 1))
            }
            lamratiominus1 <- expandSparse(lamratiominus1,
                                           n = dim(ddS)[3], across = 3)
            ##            ddSintegrand <- evalSparse3Dentrywise(lamarray * (momafter* lamratiominus1 + momchange))
            if(!hasInf) {
              #' cif is positive 
              ddSintegrand <- lamarray * (momafter* lamratiominus1 + momchange)
            } else {
              #' cif has zeroes
              zerobefore <- mapSparseEntries(ddS, 1, zerocifB,
                                             conform=TRUE, across=3)
              zerochange <- mapSparseEntries(ddS, 1, deltaInf,
                                             conform=TRUE, across=2)
              zerochange[,isdataB,] <- - zerochange[,isdataB,]
              zeroafter <- zerobefore + zerochange
              momZbefore <- mombefore
              momZafter  <- momafter
              momZbefore[ , zerocifB , ] <- 0
              IJK <- SparseIndices(zeroafter)
              momZafter[IJK] <- 0
              momZchange <- momZafter - momZbefore
              ddSintegrand <- lamarray*(momZafter*lamratiominus1 + momZchange) 
            }
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
    ## {{{{{{{{{{{{{   E N D   O F   L O O P   }}}}}}}}}}}}}

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
    geomsmooth <- geomsmooth && all(h[!isdata & ok] >= 0)
    if(mt)
      h <- data.frame(leverage=h, type=marks(loc))
    levval <- (loc %mark% h)[ok]
    levvaldum <- levval[!isdata[ok]]
    if(!mt) {
      levsmo <- Smooth(levvaldum,
                       sigma=smallsigma,
                       geometric=geomsmooth,
                       dimyx=dimyx, eps=eps)
      levnearest <- nnmark(levvaldum, dimyx=dimyx, eps=eps)
    } else {
      levsplitdum <- split(levvaldum, reduce=TRUE)
      levsmo <- Smooth(levsplitdum,
                       sigma=smallsigma,
                       geometric=geomsmooth,
                       dimyx=dimyx, eps=eps)
      levnearest <- solapply(levsplitdum, nnmark, dimyx=dimyx, eps=eps)
    }
    ## mean level
    if(fit.is.poisson) {
      a <- area(Window(loc)) * markspace.integral(loc)
      levmean <- p/a
    } else {
      levmean <- if(!mt) mean(levnearest) else mean(sapply(levnearest, mean))
    }
    lev <- list(val=levval, smo=levsmo, ave=levmean, nearest=levnearest)
    result$lev <- lev
  }
  # .......... influence .............
  if("influence" %in% what) {
    if(logi){
      X <- loc
      effX <- as.numeric(isdata) * eff - mom * (inside * logiprob)
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
      M <- as.numeric(isdata) * eff - mom * (inside * logiprob)
      M <- t(invhess %*% t(M))
      Mdum <- M
      Mdum[isdata,] <- 0
      Mdum <- Mdum / w.quad(Q)
      DFB <- msr(Q, M[isdata, ], Mdum)
    } else {
      vex <- invhess %*% t(mom)
      dex <- invhess %*% t(eff)
      switch(method,
             interpreted = {
               dis <- con <- matrix(0, nloc, ncol(mom))
               for(i in seq(nloc)) {
                 vexi <- vex[,i, drop=FALSE]
                 dexi <- dex[,i, drop=FALSE]
                 dis[i, ] <- if(isdata[i]) dexi else 0
                 con[i, ] <- if(inside[i]) (- lam[i] * vexi) else 0
               }
             },
             C = {
               dis <- t(dex)
               dis[!isdata,] <- 0
               con <- - lam * t(vex)
               con[(lam == 0 | !inside), ] <- 0
             })
      colnames(dis) <- colnames(con) <- colnames(mom)
      DFB <- msr(Q, dis[isdata, ], con)
    }
    #' add smooth component
    DFB <- augment.msr(DFB, sigma=smallsigma, dimyx=dimyx, eps=eps)
    result$dfbetas <- DFB
  }
  return(result)
}

## >>>>>>>>>>>>>>>>>>>>>>>  HELPER FUNCTIONS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

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

## >>>>>>>>>>>>>>>>  PLOT METHODS <<<<<<<<<<<<<<<<<<<<<

plot.leverage.ppm <- function(x, ...,
                              what=c("smooth", "nearest", "exact"),
                              showcut=TRUE,
                              args.cut=list(drawlabels=FALSE),
                              multiplot=TRUE) {
  what <- match.arg(what)
  fitname <- x$fitname
  defaultmain <- paste("Leverage for", fitname)

  y <- x$lev

  if(what == "exact") {
    #' plot exact quadrature locations and leverage values
    z <- do.call(plot,
                 resolve.defaults(list(x=y$val, multiplot=multiplot),
                                  list(...),
                                  list(main=defaultmain)))
    return(invisible(z))
  }

  smo <- as.im(x, what=what)
  if(is.null(smo)) return(invisible(NULL))
  
  ave <- y$ave
  if(!multiplot && inherits(smo, "imlist")) {
    ave <- ave * length(smo)
    smo <- Reduce("+", smo)
    defaultmain <- c(defaultmain, "(sum over all types of point)")
  }
  args.contour <- resolve.defaults(args.cut, list(levels=ave))
  cutinfo <- list(addcontour=showcut,
                  args.contour=args.contour)
  if(is.im(smo)) {
    do.call(plot.im,
            resolve.defaults(list(smo),
                             cutinfo,
                             list(...),
                             list(main=defaultmain)))
  } else if(inherits(smo, "imlist")) {
    do.call(plot.solist,
            resolve.defaults(list(smo),
                             cutinfo,
                             list(...),
                             list(main=defaultmain)))
  } 
  invisible(NULL)
}


persp.leverage.ppm <- function(x, ..., what=c("smooth", "nearest"),
                               main, zlab="leverage") {
  if(missing(main)) main <- deparse(substitute(x))
  what <- match.arg(what)
  y <- as.im(x, what=what)
  if(is.null(y)) return(invisible(NULL))
  if(is.im(y)) return(persp(y, main=main, ..., zlab=zlab))
  pa <- par(ask=TRUE)
  lapply(y, persp, main=main, ..., zlab=zlab)
  par(pa)
  return(invisible(NULL))
}
  
contour.leverage.ppm <- function(x, ...,
                                 what=c("smooth", "nearest"),
                                 showcut=TRUE,
                                 args.cut=list(col=3, lwd=3, drawlabels=FALSE),
                                 multiplot=TRUE) {
  defaultmain <- paste("Leverage for", x$fitname)
  smo <- as.im(x, what=what)
  y <- x$lev
  ave <- y$ave
  if(!multiplot && inherits(smo, "imlist")) {
    ave <- ave * length(smo)
    smo <- Reduce("+", smo)
    defaultmain <- c(defaultmain, "(sum over all types of point)")
  }
  
  argh1 <- resolve.defaults(list(...),
                            list(main=defaultmain))
  argh2 <- resolve.defaults(args.cut,
                            list(levels=ave),
                            list(...))

  if(is.im(smo)) {
    #' single panel
    out <- do.call(contour, append(list(x=smo), argh1))
    if(showcut)
      do.call(contour, append(list(x=smo, add=TRUE), argh2))
  } else if(inherits(smo, "imlist")) {
    #' multiple panels
    argh <- append(list(x=smo, plotcommand ="contour"), argh1)
    if(showcut) {
      argh <- append(argh,
                     list(panel.end=function(i, y, ...) contour(y, ...),
                          panel.end.args=argh2))
    } 
    out <- do.call(plot.solist, argh) 
  } else {
    warning("Unrecognised format")
    out <- NULL
  }
  return(invisible(out))
}

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

## >>>>>>>>>>>>>>>>  CONVERSION METHODS <<<<<<<<<<<<<<<<<<<<<


as.im.leverage.ppm <- function(X, ..., what=c("smooth", "nearest")) {
  what <- match.arg(what)
  y <- switch(what,
              smooth = X$lev$smo,
              nearest = X$lev$nearest)
  if(is.null(y))
    warning(paste("Data for", sQuote(what), "image are not available:",
                  "please recompute the leverage using the current spatstat"),
            call.=FALSE)
  return(y) # could be either an image or a list of images
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


## >>>>>>>>>>>>>>>>  PRINT METHODS <<<<<<<<<<<<<<<<<<<<<

print.leverage.ppm <- function(x, ...) {
  splat("Point process leverage function")
  fitname <- x$fitname
  splat("for model:", fitname)
  lev <- x$lev
  splat("\nExact values:")
  print(lev$val)
  splat("\nSmoothed values:")
  print(lev$smo)
  ## for compatibility we retain the x$fit usage
  if(x$fit.is.poisson %orifnull% is.poisson(x$fit))
    splat("\nAverage value:", lev$ave)
  return(invisible(NULL))
}

print.influence.ppm <- function(x, ...) {
  splat("Point process influence measure")  
  fitname <- x$fitname
  splat("for model:", fitname)
  splat("\nExact values:")
  print(x$infl)
  return(invisible(NULL))
}

## >>>>>>>>>>>>>>>>  SUBSET METHODS <<<<<<<<<<<<<<<<<<<<<

"[.leverage.ppm" <- function(x, i, ..., update=TRUE) {
  if(missing(i)) return(x)
  y <- x$lev
  smoi <- if(is.im(y$smo)) y$smo[i, ...] else solapply(y$smo, "[", i=i, ...)
  if(!inherits(smoi, c("im", "imlist"))) return(smoi)
  y$smo <- smoi
  y$val <- y$val[i, ...]
  if(update) 
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

## >>>>>>>>>>>>>>>>  GEOMETRICAL OPERATIONS <<<<<<<<<<<<<<<<<<<<<

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

