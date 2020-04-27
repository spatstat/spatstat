#
#    relrisk.lpp.R
#
#   Estimation of relative risk on network
#
#  $Revision: 1.6 $  $Date: 2020/04/27 03:08:26 $
#

relrisk.lpp <- local({

  relrisk.lpp <- function(X, sigma, ..., 
                          at=c("pixels", "points"), 
                          relative=FALSE,
                          adjust=1, 
                          casecontrol=TRUE, control=1, case,
                          finespacing=FALSE) {
    stopifnot(is.lpp(X))
    stopifnot(is.multitype(X))
    control.given <- !missing(control)
    case.given <- !missing(case)
    at <- match.arg(at)
    ##
    marx <- marks(X)
    types <- levels(marx)
    ntypes <- length(types)
    ## 
    if(ntypes == 1L)
      stop("Data contains only one type of points")
    casecontrol <- casecontrol && (ntypes == 2L)
    if((control.given || case.given) && !(casecontrol || relative)) {
      aa <- c("control", "case")[c(control.given, case.given)]
      nn <- length(aa)
      warning(paste(ngettext(nn, "Argument", "Arguments"),
                    paste(sQuote(aa), collapse=" and "),
                    ngettext(nn, "was", "were"),
                    "ignored, because relative=FALSE and",
                    if(ntypes==2L) "casecontrol=FALSE" else
                    "there are more than 2 types of points"))
    }
    ## compute bandwidth
    if(is.function(sigma)) {
      sigma <- do.call.matched(sigma, list(X=X, ...))
      if(!is.numeric(sigma))
        stop("The function 'sigma' did not return a numerical value",
             call.=FALSE)
    }
    check.1.real(sigma) # includes Inf
    sigma <- adjust * as.numeric(sigma)
    ## .........................................
    ## compute intensity estimates for each type
    ## .........................................
    Y <- split(X)
    switch(at,
           pixels = {
             ## intensity estimates of each type
             Deach <- solapply(Y, density.lpp, sigma=sigma,
                               ..., finespacing=finespacing)
             ## compute intensity estimate for unmarked pattern
             Dall  <- density(unmark(X), sigma=sigma,
                              ..., finespacing=finespacing)
           },
           points = {
             ## intensity estimates of each type **at each data point**
             Deachfun <- solapply(Y, densityfun.lpp, sigma=sigma,
                                  ..., finespacing=finespacing)
             Deach <- as.data.frame(sapply(Deachfun, function(f, P) f(P), P=X))
             ## leave-one-out estimates
             Dself <- lapply(Y, density.lpp, sigma=sigma,
                             at="points", leaveoneout=TRUE,
                             ..., finespacing=finespacing)
             ## insert leave-one-out estimates in correct place
             Deachsplit <- split(Deach, marx)
             for(j in 1:ntypes) {
               Deachsplit[[j]][, j] <- Dself[[j]]
             }
             split(Deach, marx) <- Deachsplit
             ## total
             Dall <- rowSums(Deach)
           })
    ## .........................................
    ## compute probabilities/risks
    ## .........................................
    if(ntypes == 2 && casecontrol) {
      if(control.given || !case.given) {
        stopifnot(length(control) == 1)
        if(is.numeric(control)) {
          icontrol <- control <- as.integer(control)
          stopifnot(control %in% 1:2)
        } else if(is.character(control)) {
          icontrol <- match(control, types)
          if(is.na(icontrol)) stop(paste("No points have mark =", control))
        } else
          stop(paste("Unrecognised format for argument", sQuote("control")))
        if(!case.given)
          icase <- 3 - icontrol
      }
      if(case.given) {
        stopifnot(length(case) == 1)
        if(is.numeric(case)) {
          icase <- case <- as.integer(case)
          stopifnot(case %in% 1:2)
        } else if(is.character(case)) {
          icase <- match(case, types)
          if(is.na(icase)) stop(paste("No points have mark =", case))
        } else stop(paste("Unrecognised format for argument", sQuote("case")))
        if(!control.given) 
          icontrol <- 3 - icase
      }
      ## compute ......
      switch(at,
             pixels = {
               ## compute probability of case
               pcase <- Deach[[icase]]/Dall
               ## correct small numerical errors
               pcase <- clamp01(pcase)
               ## trap NaN values
               nbg <- badvalues(pcase)
               if(any(nbg)) {
                 ## apply l'Hopital's rule:
                 ##     p(case) = 1{nearest neighbour is case}
                 distcase <- as.linim(distfun(Y[[icase]]))
                 distcontrol <- as.linim(distfun(Y[[icontrol]]))
                 closecase <- eval.linim(as.integer(distcase < distcontrol))
                 pcase[nbg] <- closecase[nbg]
               }
               if(!relative) {
                 result <- pcase
               } else {
                 result <- eval.im(ifelse(pcase < 1, pcase/(1-pcase), NA))
               }
             },
             points={
               ## compute probability of case
               pcase <- Deach[,icase]/Dall
               ## correct small numerical errors
               pcase <- clamp01(pcase)
               ## trap NaN values
               if(any(nbg <- badvalues(pcase))) {
                 ## apply l'Hopital's rule
                 imarks <- as.integer(marx)
                 nntype <- imarks[nnwhich(X)]
                 pcase[nbg] <- as.integer(nntype[nbg] == icase)
               }
               if(!relative) {
                 result <- pcase
               } else {
                 result <- ifelse(pcase < 1, pcase/(1-pcase), NA)
               }
             })
    } else {
      ## several types
      if(relative) {
        ## need 'control' type
        stopifnot(length(control) == 1)
        if(is.numeric(control)) {
          icontrol <- control <- as.integer(control)
          stopifnot(control %in% 1:ntypes)
        } else if(is.character(control)) {
          icontrol <- match(control, types)
          if(is.na(icontrol)) stop(paste("No points have mark =", control))
        } else
          stop(paste("Unrecognised format for argument", sQuote("control")))
      }
      switch(at,
             pixels={
               probs <- as.solist(lapply(Deach, "/", e2=Dall))
               ## correct small numerical errors
               probs <- as.solist(lapply(probs, clamp01))
               ## trap NaN values
               nbg <- lapply(probs, badvalues)
               nbg <- Reduce("|", nbg)
               if(any(nbg)) {
                 ## apply l'Hopital's rule
                 whichnn <- as.linim(nnfun(X))
                 imarks <- as.integer(marx)
                 typenn <- eval.im(imarks[whichnn])
                 typennsub <- typenn[nbg]
                 for(k in seq_along(probs)) 
                   probs[[k]][nbg] <- (typennsub == k)
               }
               if(!relative) {
                 result <- probs
               } else {
                 result <- solapply(probs,
                                    divideifpositive,
                                    d = probs[[icontrol]])
               }
             },
             points = {
               probs <- Deach/Dall
               ## correct small numerical errors
               probs <- clamp01(probs)
               ## trap NaN values
               bad <- badvalues(probs)
               badrow <- matrowany(bad)
               if(any(badrow)) {
                 ## apply l'Hopital's rule
                 imarks <- as.integer(marx)
                 typenn <- imarks[nnwhich(X)]
                 probs[badrow, ] <- (typenn == col(result))[badrow, ]
               }
               if(!relative) {
                 result <- probs
               } else {
                 result <- probs/probs[,icontrol]
               }
            })
    }
    attr(result, "sigma") <- sigma
    return(result)
  }

  clamp01 <- function(x) {
    if(is.linim(x)) return(eval.linim(pmin(pmax(x, 0), 1)))
    if(is.im(x)) return(eval.im(pmin(pmax(x, 0), 1)))
    if(is.data.frame(x)) x <- as.matrix(x)
    return(pmin(pmax(x, 0), 1))
  }

  badvalues <- function(x) {
    if(is.linim(x)) return(eval.linim(!is.finite(x)))
    if(is.im(x)) return(eval.im(!is.finite(x)))
    if(is.data.frame(x)) x <- as.matrix(x)
    return(!(is.finite(x) | is.na(x)))
  }

  divideifpositive <- function(z, d) { eval.linim(ifelse(d > 0, z/d, NA)) }
  
  relrisk.lpp
})


bw.relrisklpp <- local({

hargnames <- c("hmin", "hmax")

bw.relrisklpp <- function(X, ..., 
                          method=c("likelihood",
                                   "leastsquares",
                                   "KelsallDiggle",
                                   "McSwiggan"),
                          distance=c("path", "euclidean"),
                          hmin=NULL, hmax=NULL,
                          nh=NULL,
                          fast=TRUE, fastmethod="onestep", floored=TRUE,
                          reference=c("thumb", "uniform", "sigma"),
                          allow.infinite=TRUE,
                          epsilon=1e-20, fudge=0,
                          verbose=FALSE, warn=TRUE) {
  startTime <- proc.time()
  stopifnot(is.lpp(X))
  stopifnot(is.multitype(X))
  method <- match.arg(method)
  reference <- match.arg(reference)
  distance <- match.arg(distance)
  if(is.null(nh)) nh <- switch(distance, path=256, euclidean=16)
  ## validate X
  marx <- marks(X)
  types <- levels(marx)
  ntypes <- length(types)
  if(ntypes == 1L)
    stop("There is only one type of point", call.=FALSE)
  if(ntypes > 2L && distance == "path")
    stop(paste("Sorry, bw.relrisklpp(distance='path') is not yet supported",
               "for > 2 types of points"), call.=FALSE)
  ## determine range of bandwidths
  if(got.hmax <- !missing(hmax)) { check.1.real(hmax) ; stopifnot(hmax > 0) }
  if(got.hmin <- !missing(hmin)) { check.1.real(hmin) ; stopifnot(hmin > 0) }
  if(got.hmax && got.hmin) {
    stopifnot(hmin < hmax)
  } else if(got.hmax) {
    hmin <- hmax/20
  } else if(got.hmin) {
    hmax <- hmin * 20
  } else {
    ss <- bw.scott.iso(X)
    dd <- diameter(Frame(X))                                               
    srange <- range(c(ss/10, ss*5, dd/5))
    hmin <- srange[1L]
    hmax <- srange[2L]
  }
  if(verbose) splat("Bandwidth range:", prange(c(hmin, hmax)))
  ##
  if(distance == "euclidean") {
    if(verbose) splat("Euclidean smoothing")
    if(method %in% c("McSwiggan", "KelsallDiggle"))
      stop(paste0("Sorry, bw.relrisklpp(method=", sQuote(method),
                  ") is not yet supported for > 2 types of points"),
           call.=FALSE)
    sigmavalues <- seq(hmin, hmax, length.out=nh)
    cv <- numeric(nh)
    witch <- cbind(seq_along(marx), as.integer(marx))
    pstate <- list()
    if(verbose) cat(paste("Processing", nh, "values of bandwidth ..."))
    for(i in 1:nh) {
      si <- sigmavalues[i]
      p <- relrisk(X, si, at="points", distance="euclidean", casecontrol=FALSE)
      pobs <- p[witch]
      cv[i] <- switch(method,
                      likelihood=log(prod(pobs)),
                      leastsquares=sum((1-pobs)^2))
      if(verbose) pstate <- progressreport(i, nh, state=pstate)
    }
    result <- switch(method,
                     likelihood = bw.optim(cv, sigmavalues, optimum="max", 
                              hname="sigma", cvname="logL",  
                              criterion="likelihood cross-validation",
                              hargnames=hargnames,
                              unitname=unitname(X)),
                     leastsquares = bw.optim(cv, sigmavalues, 
                              hname="sigma", cvname="psq", 
                              criterion="least squares cross-validation",
                              hargnames=hargnames,
                              unitname=unitname(X)))
    return(result)
  }
  ## ---------- heat kernel (distance='path') ------------------------------
  sigma <- hmax
  nsigma <- ceiling(nh * hmax/(hmax-hmin))
  #'
  if(verbose) splat("Setting up network data...")
  L <- domain(X)
  TOTLEN <- volume(L)
  g <- densityfun.lpp(X=unmark(X), sigma=sigma, nsigma=nsigma,
                      exit="setup", verbose=FALSE, ...)
  #' extract internal data
  finenet     <- g$linnet_obj
  lixelmap    <- g$lixelmap
  lixelweight <- g$lixelweight
  Amatrix     <- g$Amatrix
  ## U0          <- g$U0  # not used
  deltax      <- g$deltax
  deltat      <- g$deltat
  #'
  if(allow.infinite) {
    df <- as.data.frame(vertices(finenet))[,c("x","y","segcoarse","tpcoarse")]
    colnames(df) <- c("x", "y", "seg", "tp")
    fineverticescoarsenet <- lpp(df, L)
  }
  ## split into types
  Y <- split(X)
  X1 <- Y[[1L]]
  X2 <- Y[[2L]]
  n1 <- npoints(X1)
  n2 <- npoints(X2)
  #' discretise X1, X2 separately
  #' Each data point is mapped to two endpoints of a tiny segment
  I1 <- (as.integer(marx) == 1L)
  lixelweight1 <- lixelweight[I1]
  lixelmap1     <- lixelmap[I1]
  U01 <- tapplysum(lixelweight1, list(lixelmap1))
  I2 <- !I1
  lixelweight2 <- lixelweight[I2]
  lixelmap2     <- lixelmap[I2]
  U02 <- tapplysum(lixelweight2, list(lixelmap2))
  #' determine number of time steps
  niter <- round((sigma^2)/(2 * deltat))
  nsample <- length(U01)
  #' solve heat equation separately for X1 and X2
  if(verbose) splat("Computing intensity estimates",
                    "for", nsigma, "out of",
                    niter, "bandwidth values at",
                    nsample, "sample locations ...")
  K1 <- K2 <- matrix(0, nsample, nsigma)
  U1 <- U01
  U2 <- U02
  blocksize <- ceiling(niter/nsigma)
  pstate <- list()
  for(isave in 1:nsigma) {
    nit <- min(blocksize, niter - (isave-1L)*blocksize)
    if(nit > 0) {
      for(iter in 1:nit) {
        U1 <- as.numeric(Amatrix %*% U1)
        U2 <- as.numeric(Amatrix %*% U2)
      }
      if(verbose) pstate <- progressreport(isave, nsigma, state=pstate)
    }
    K1[, isave] <- U1 
    K2[, isave] <- U2
  }
  if(verbose) splat("Done.")
  #' add small amount to log intensity
  logK1 <- log(K1 + epsilon)
  logK2 <- log(K2 + epsilon)
  logK1 <- t(logK1)
  logK2 <- t(logK2)

  #' Map each data point to closest endpoint
  J1 <- closeroftwo(lixelweight1, lixelmap1)
  J2 <- closeroftwo(lixelweight2, lixelmap2)
  
  #' Term ghat from Term 2 - intensity of Type 2 events at Type 1 locations
  ghat <- K2[J1, ] + epsilon
  ghat <- t(ghat)

  #' Term fhat from Term 3 - intensity of Type 1 events at Type 2 locations
  fhat <- K1[J2,] + epsilon
  fhat <- t(fhat)

  #' For Term 2 calculate f^(-i)_{h_1}(x_i)
  #' = intensity at type 1 event x_i estimated from all type 1 events except x_i
  #' Likewise g^(-j)_{h_1}(y_j)
  #' = intensity at type 2 event y_j estimated from all type 2 events except y_j
  if(verbose) splat("Computing leave-one-out estimates at data points.")
  if(verbose) cat("Type 1 ...")
  fminusi <- densitypointsLPP(X1, sigma, dx=deltax, dt=deltat,
                              nsigma=nsigma, leaveoneout=TRUE,
                              fast=fast, fastmethod=fastmethod,
                              floored=floored)
  if(verbose) cat(" Done.\nType 2 ...")
  gminusj <- densitypointsLPP(X2, sigma, dx=deltax, dt=deltat,
                              nsigma=nsigma, leaveoneout=TRUE,
                              fast=fast, fastmethod=fastmethod,
                              floored=floored)
  fminusi <- t(fminusi)
  gminusj <- t(gminusj)
  tau <- attr(gminusj, "sigma")
  use <- (hmin <= tau) & (tau <= hmax)
  if(verbose) splat("Done.")

  #' reference intensity (used in McSwiggan (modified K-D) method)
  switch(reference,
         sigma = {
           #' Use largest value of sigma
           #' reference intensity of type 1 process
           fbar <- K1[,nsigma] + epsilon
           #' reference intensity of type 2 process
           gbar <- K2[, nsigma] + epsilon
           #' leave-one-out estimates at data points
           fbarminusi <- fminusi[nsigma, ]
           gbarminusj <- gminusj[nsigma, ]
         },
         thumb = {
           #' Use smoothers selected by rule of thumb
           b1 <- bw.scott.iso(X1)
           b2 <- bw.scott.iso(X2)
           i1 <- which.min(abs(b1 - tau))
           i2 <- which.min(abs(b2 - tau))
           #' reference intensity of type 1 process           
           fbar <- K1[,i1] + epsilon
           #' reference intensity of type 2 process           
           gbar <- K2[,i2] + epsilon
           #' leave-one-out estimates at data points
           fbarminusi <- fminusi[i1, ]
           gbarminusj <- gminusj[i2, ]
         },
         uniform = {
           #' Use uniform intensity
           fbar <- rep.int(n1/TOTLEN, nrow(K1))
           gbar <- rep.int(n2/TOTLEN, nrow(K2))
           #' leave-one-out estimates at data points
           fbarminusi <- rep.int((n1-1)/TOTLEN, n1)
           gbarminusj <- rep.int((n2-1)/TOTLEN, n2)
         })
  #'   reference intensity of type 1 process at type 1 points
  #' fbari <- fbar[J1]  # not used
  #'   reference intensity of type 2 process at type 2 points
  #' gbarj <- gbar[J2]  # not used
  
  #' Avoid very small estimates
  if(fudge > 0) {
    minloo <- fudge/TOTLEN
    gminusj[] <- pmax(minloo, gminusj[])
    fminusi[] <- pmax(minloo, fminusi[])
    gbarminusj[] <- pmax(minloo, gbarminusj[])
    fbarminusi[] <- pmax(minloo, fbarminusi[])
  } else {
    gminusj <- gminusj + epsilon
    fminusi <- fminusi + epsilon
    gbarminusj <- gbarminusj + epsilon
    fbarminusi <- fbarminusi + epsilon
  }

  if(allow.infinite) {
    ## also compute values for sigma = Inf
    ## corresponding to a constant relative risk
    if(verbose) splat("Computing estimates for sigma=Inf ...")
    fInfFun <- densityfun(X1, Inf)
    gInfFun <- densityfun(X2, Inf)
    fhatInf <- fInfFun(X2) # intensity of X1 at points of X2
    ghatInf <- gInfFun(X1) # intensity of X2 at points of X1
    K1Inf <- fInfFun(fineverticescoarsenet) # intensity of X1 at fine grid
    K2Inf <- gInfFun(fineverticescoarsenet) # intensity of X2 at fine grid
    #' intensity of X1 at points of X1, leave-one-out
    fminusiInf <- densitypointsLPP(X1, Inf, leaveoneout=TRUE)
    gminusjInf <- densitypointsLPP(X2, Inf, leaveoneout=TRUE)
    #' ensure they are row vectors 
    fhatInf <- matrix(fhatInf[], nrow=1)
    ghatInf <- matrix(ghatInf[], nrow=1)
    fminusiInf <- matrix(fminusiInf[], nrow=1)
    gminusjInf <- matrix(gminusjInf[], nrow=1)
    K1Inf <- matrix(K1Inf[], nrow=1)
    K2Inf <- matrix(K2Inf[], nrow=1)
    logK1Inf <- log(K1Inf + epsilon)
    logK2Inf <- log(K2Inf + epsilon)
  }
  
  #' Compute terms in cross-validation score
  if(verbose) splat("Computing basic cross-validation terms ...")
  Term1 <- deltax * xvalterm1(logK1,logK2)
  Term1Inf <- deltax * xvalterm1(logK1Inf,logK2Inf)

  switch(method,
         KelsallDiggle = {
           #' ........... original Kelsall-Diggle criterion .................
           if(verbose) splat("Computing Kelsall-Diggle criterion ...")
           Term2 <- (-2) * xvalterm2(fminusi, ghat) 
           Term3 <- (-2) * xvalterm2(gminusj, fhat)
           ## Term3 <- t(Term3)
           CKD <- -Term1 + Term2 + Term3
           #' repeat for sigma=Inf
           Term2Inf <- (-2) * xvalterm2(fminusiInf, ghatInf) 
           Term3Inf <- (-2) * xvalterm2(gminusjInf, fhatInf)
           ## Term3Inf <- t(Term3Inf)
           Cinf <- -Term1Inf + Term2Inf + Term3Inf
           ##
           CKDout <- c(CKD[use], Cinf)
           tauout <- c(tau[use], Inf)
           result <- bw.optim(CKDout, tauout,
                              hname="sigma", cvname="C",
                              criterion="Kelsall-Diggle cross-validation",
                              hargnames=hargnames,
                              unitname=unitname(X))
         },
         McSwiggan = {
           #' .............. modified criterion .....................
           if(verbose) splat("Computing modified Kelsall-Diggle criterion ...")
           ModTerm2 <- -2 * xvalterm4(fminusi, ghat, 1/fbarminusi)
           ModTerm3 <- -2 * xvalterm4(gminusj, fhat, 1/gbarminusj)
           Term4 <- -2 * deltax * xvalterm4(t(K1), t(K2), log(fbar/gbar))
           Term4[!is.finite(Term4)] <- Inf
           modC <- Term1 + ModTerm2 + ModTerm3 + Term4
           ## again for sigma=Inf
           ModTerm2Inf <- -2 * xvalterm4(fminusiInf, ghatInf, 1/fbarminusi)
           ModTerm3Inf <- -2 * xvalterm4(gminusjInf, fhatInf, 1/gbarminusj)
           Term4Inf <- -2 * deltax * xvalterm4(K1Inf, K2Inf, log(fbar/gbar))
           Term4Inf[!is.finite(Term4Inf)] <- Inf
           Cinf <- Term1Inf + ModTerm2Inf + ModTerm3Inf + Term4Inf
           modCout <- c(modC[use], Cinf)
           tauout <- c(tau[use], Inf)
           result <- bw.optim(modCout, tauout,
                              hname="sigma", cvname="Cmod",
                criterion="McSwiggan modified Kelsall-Diggle cross-validation",
                hargnames=hargnames,
                unitname=unitname(X))
         },
         likelihood = {
           #' .............. likelihood criterion .....................
           if(verbose) splat("Computing likelihood criterion ...")
           TermA <- xvalterm5(fminusi, ghat)
           TermB <- xvalterm5(gminusj, fhat)
           loglik <- TermA + TermB
           ##
           TermAInf <- xvalterm5(fminusiInf, ghatInf)
           TermBInf <- xvalterm5(gminusjInf, fhatInf)
           loglikInf <- TermAInf + TermBInf
           ##
           loglikout <- c(loglik[use], loglikInf)
           tauout    <- c(tau[use], Inf)
           # as.numeric(loglikInf)
           result <- bw.optim(loglikout, tauout, optimum="max", 
                              hname="sigma", cvname="logL",  
                              criterion="likelihood cross-validation",
                              hargnames=hargnames,
                              unitname=unitname(X))
         },
         leastsquares = {
           #' .............. least squares criterion .....................
           if(verbose) splat("Computing least squares criterion ...")
           Term6A <- xvalterm6(fminusi, ghat)
           Term6B <- xvalterm6(gminusj, fhat)
           sqprob <- Term6A + Term6B
           #' 
           Term6AInf <- xvalterm6(fminusiInf, ghatInf)
           Term6BInf <- xvalterm6(gminusjInf, fhatInf)
           sqprobInf <- Term6AInf + Term6BInf
           ##
           sqprobout <- c(sqprob[use], sqprobInf)
           tauout    <- c(tau[use], Inf)
           result <- bw.optim(sqprobout, tauout, 
                              hname="sigma", cvname="psq", 
                              criterion="least squares cross-validation",
                              hargnames=hargnames,
                              unitname=unitname(X))
         })
  if(verbose) splat("Done.")
  result <- timed(result, starttime=startTime)
  return(result)
}  

closeroftwo <- function(ww, ff) {
  even <- c(FALSE,TRUE)
  odd  <- c(TRUE, FALSE)
  as.integer(ifelse(ww[even] > ww[odd], ff[even], ff[odd]))
}

xvalterm1 <- function(x, y) { rowSums((x-y)^2) }
xvalterm2 <- function(x, y) { rowSums((log(x/y))/x) }
xvalterm4 <- function(x, y, w) { as.numeric(log(x/y) %*% w) }
xvalterm5 <- function(x, y) { rowSums(log(x/(x+y))) }
xvalterm6 <- function(x, y) { rowSums((1 - x/(x+y))^2) }

bw.relrisklpp


})



                          
