##
## Asymptotic covariance & correlation matrices
## and Fisher information matrix
## for ppm objects
##
##  $Revision: 1.111 $  $Date: 2014/10/24 00:22:30 $
##

vcov.ppm <- local({

vcov.ppm <- function(object, ..., what="vcov", verbose=TRUE,
                     gam.action=c("warn", "fatal", "silent"),
                     matrix.action=c("warn", "fatal", "silent"),
                     logi.action=c("warn", "fatal", "silent"),
                     hessian=FALSE) {
  verifyclass(object, "ppm")
  argh <- list(...)

  gam.action <- match.arg(gam.action)
  matrix.action <- match.arg(matrix.action)
  logi.action <- match.arg(logi.action)

  stopifnot(length(what) == 1 && is.character(what))
  what.options <- c("vcov", "corr", "fisher", "Fisher", "internals", "all")
  what.map     <- c("vcov", "corr", "fisher", "fisher", "internals", "all")
  if(is.na(m <- pmatch(what, what.options)))
    stop(paste("Unrecognised option: what=", sQuote(what)))
  what <- what.map[m]

  ## no parameters, no variance
  if(length(coef(object)) == 0) {
    result <- switch(what,
                     vcov=, corr=, fisher= {
                       matrix(, 0, 0)
                     },
                     internals=, all={
                       list()
                     })
    return(result)
  }
  
  ## nonstandard calculations (hack) 
  generic.triggers <- c("A1", "A1dummy", "new.coef", "matwt", "saveterms")
  nonstandard <- any(generic.triggers %in% names(argh))
  saveterms <- identical(resolve.1.default("saveterms", argh), TRUE)
  
  ## Fisher information *may* be contained in object
  fisher <- object$fisher
  varcov <- object$varcov
  
  ## Do we need to go into the guts?
  needguts <- nonstandard ||
    (is.null(fisher) && what=="fisher") ||
    (is.null(varcov) && what %in% c("vcov", "corr")) ||
    (what %in% c("internals", "all")) 

  ## In general it is not true that varcov = solve(fisher)
  ## because we might use different estimators,
  ## or the parameters might be a subset of the canonical parameter

  if(needguts) {
    ## warn if fitted model was obtained using GAM
    if(identical(object$fitter, "gam")) {
      switch(gam.action,
             fatal={
               stop(paste("model was fitted by gam();",
                          "execution halted because fatal=TRUE"),
                    call.=FALSE)
             },
             warn={
               warning(paste("model was fitted by gam();",
                             "asymptotic variance calculation ignores this"),
                       call.=FALSE)
             },
             silent={})
    }
    ## ++++ perform main calculation ++++
    if((is.poisson(object) || hessian) && object$method != "logi") {
      ## Poisson model, or Hessian of Gibbs model
      results <- vcalcPois(object, ..., what=what,
                           matrix.action=matrix.action,
                           verbose=verbose, fisher=fisher)
    } else {
      ## Gibbs model 
      results <- vcalcGibbs(object, ..., what=what,
                            matrix.action=matrix.action)
    }
    varcov <- results$varcov
    fisher <- results$fisher
    internals  <- results$internals
  }
  
  if(what %in% c("vcov", "corr") && is.null(varcov)) {
    ## Need variance-covariance matrix.
    if(!is.null(fisher) && is.poisson(object)) 
      ## Derive from Fisher information
      varcov <- checksolve(fisher, matrix.action,
                           "Fisher information matrix",
                           "variance")
  }

  out <- switch(what,
                fisher = fisher,
                vcov   = varcov,
                corr   = {
                  if(is.null(varcov)) return(NULL)
                  sd <- sqrt(diag(varcov))
                  varcov / outer(sd, sd, "*")
                },
                internals = internals,
                all = results
                )
  return(out)
}

## ................  variance calculation for Poisson models  .............

vcalcPois <- function(object, ...,
                      what = c("vcov", "corr", "fisher", "internals", "all"),
                      matrix.action=c("warn", "fatal", "silent"),
                      method=c("C", "interpreted"),
                      verbose=TRUE,
                      fisher=NULL, 
                      matwt=NULL, new.coef=NULL,
                      saveterms=FALSE) {
  ## variance-covariance matrix of Poisson model,
  ## or Hessian of Gibbs model
  what <- match.arg(what)
  method <- match.arg(method)
  matrix.action <- match.arg(matrix.action)
  if(reweighting <- !is.null(matwt)) 
    stopifnot(is.numeric(matwt) && is.vector(matwt))
  internals <- NULL
  nonstandard <- reweighting || !is.null(new.coef) || saveterms
  ## compute Fisher information if not known
  if(is.null(fisher) || nonstandard) {
    gf <- getglmfit(object)
    ## we need a glm or gam
    if(is.null(gf)) {
      if(verbose) 
        warning("Refitting the model using GLM/GAM")
      object <- update(object, forcefit=TRUE)
      gf <- getglmfit(object)
      if(is.null(gf))
        stop("Internal error - refitting did not yield a glm object")
    }
    ## compute fitted intensity and sufficient statistic
    ltype <- if(is.poisson(object)) "trend" else "lambda"
    lambda <- fitted(object, type=ltype, new.coef=new.coef, check=FALSE)
    mom <- model.matrix(object)
    nmom <- nrow(mom)
    Q <- quad.ppm(object)
    wt <- w.quad(Q)
    ok <- getglmsubset(object)
    Z  <- is.data(Q)
    ## save them
    if(what == "internals") {
      internals <-
        if(!saveterms) list(suff=mom) else
      list(suff=mom, mom=mom, lambda=lambda, Z=Z, ok=ok)
    }
    ## Now restrict all terms to the domain of the pseudolikelihood
    lambda <- lambda[ok]
    mom <- mom[ok, , drop=FALSE]
    wt <- wt[ok]
    Z <- Z[ok]
    ## apply weights to rows of model matrix - temporary hack
    if(reweighting) {
      nwt <- length(matwt)
      if(nwt == nmom) {
        ## matwt matches original quadrature scheme - trim it
        matwt <- matwt[ok]
      } else if(nwt != sum(ok))
        stop("Hack argument matwt has incompatible length")
      mom.orig <- mom
      mom <- matwt * mom
    }
    ## compute Fisher information
    switch(method,
           C = {
             fisher <- sumouter(mom, lambda * wt)
             if(reweighting) {
               gradient <- sumouter(mom.orig, matwt * lambda * wt)
             }
           },
           interpreted = {
             if(!reweighting) {
               fisher <- 0
               for(i in 1:nrow(mom)) {
                 ro <- mom[i, ]
                 v <- outer(ro, ro, "*") * lambda[i] * wt[i]
                 if(!any(is.na(v)))
                   fisher <- fisher + v
               }
               momnames <- dimnames(mom)[[2]]
               dimnames(fisher) <- list(momnames, momnames)
             } else {
               fisher <- gradient <- 0
               for(i in 1:nrow(mom)) {
                 ro <- mom[i, ]
                 ro0 <- mom.orig[i,]
                 ldu <- lambda[i] * wt[i]
                 v <- outer(ro, ro, "*") * ldu
                 v0 <- outer(ro0, ro0, "*") * matwt[i] * ldu
                 if(!any(is.na(v)))
                   fisher <- fisher + v
                 if(!any(is.na(v0)))
                   gradient <- gradient + v0
               }
               momnames <- dimnames(mom)[[2]]
               dn <- list(momnames, momnames)
               dimnames(fisher) <- dimnames(gradient) <- dn
             }
           })
  } 

  if(what %in% c("all", "internals")) {
    ## Internals needed
    if(is.null(internals))
      internals <- list(suff = model.matrix(object))
    internals$fisher <- fisher
    if(reweighting)
      internals$gradient <- gradient
    ilist <- list(internals=internals)
  }

  if(what %in% c("all", "vcov", "corr")) {
    ## Variance-covariance matrix needed
    if(!reweighting) {
      ## Derive variance-covariance from Fisher info
      varcov <- checksolve(fisher, matrix.action,
                           "Fisher information matrix",
                           "variance")
      vcovlist <- list(fisher=fisher, varcov=varcov)
    } else {
      invgrad <- checksolve(gradient, matrix.action,
                            "gradient matrix", "variance")
      varcov <- if(is.null(invgrad)) NULL else
      invgrad %*% fisher %*% invgrad
      vcovlist <- list(fisher=fisher, varcov=varcov, invgrad=invgrad)
    }
  }
  result <- switch(what,
                   fisher    = list(fisher=fisher),
                   vcov      = vcovlist,
                   corr      = vcovlist,
                   internals = ilist,
                   all       = append(ilist, vcovlist))
  return(result)
}


## ...................... vcov calculation for Gibbs models ....................

vcalcGibbs <- function(fit, ...,
                       what = c("vcov", "corr", "fisher", "internals", "all"),
                       generic=FALSE) {
  what <- match.arg(what)

  if(missing(generic)) {
    ## Change default to TRUE in certain cases
    ## For logistic fits, use generic method by default
    if(fit$method == "logi")
      generic <- TRUE
    ## For 'difficult' interactions, use generic method by default
    fasterbygeneric <- c("Areainter")
    if(as.interact(fit)$creator %in% fasterbygeneric)
      generic <- TRUE
  }
  
  ## decide whether to use the generic algorithm
  generic.triggers <- c("A1", "A1dummy", "hessian",
                        "new.coef", "matwt", "saveterms")
  
  use.generic <-
    generic ||
  !is.stationary(fit) ||
  (fit$method == "logi" && ("marks" %in% variablesinformula(fit$trend))) ||
  (fit$method != "logi" && has.offset(fit)) ||
  (fit$method == "logi" && has.offset.term(fit)) ||
  !(fit$correction == "border" && fit$rbord == reach(fit)) ||
  any(generic.triggers %in% names(list(...))) ||
  !identical(options("contrasts")[[1]],
             c(unordered="contr.treatment",
               ordered="contr.poly"))
  
  ## compute
  spill <- (what %in% c("all", "internals", "fisher"))
  spill.vc <- (what == "all")
  out <- if(use.generic)
    vcalcGibbsGeneral(fit, ..., spill=spill, spill.vc=spill.vc) else
    vcalcGibbsSpecial(fit, ..., spill=spill, spill.vc=spill.vc)

  switch(what,
         vcov = ,
         corr = {
           ## out is the variance-covariance matrix; return it
           return(list(varcov=out))
         },
         fisher = {
           ## out is a list of internal data: extract the Fisher info
           Fmat <- with(out,
                        if(fit$method != "logi") Sigma else Sigma1log+Sigma2log)
           return(list(fisher=Fmat))
         },
         internals = {
           ## out is a list of internal data: return it
           ## (ensure model matrix is included)
           if(is.null(out$mom))
             out$mom <- model.matrix(fit)
           return(list(internals=out))
         },
         all = {
           ## out is a list(internals, vc): return it
           ## (ensure model matrix is included)
           if(is.null(out$internals$mom))
             out$internals$mom <- model.matrix(fit)
           ## ensure Fisher info is included
           if(is.null(out$internals$fisher)) {
             Fmat <- with(out$internals,
                     if(fit$method != "logi") Sigma else Sigma1log+Sigma2log)
             out$internals$fisher <- Fmat
           }
           return(out)
         },
         )
  return(NULL)
}

## ...................... general algorithm ...........................

vcalcGibbsGeneral <- function(model,
                         ...,
                         spill = FALSE,
                         spill.vc = FALSE,
                         matrix.action=c("warn", "fatal", "silent"),
                         logi.action=c("warn", "fatal", "silent"),
                         algorithm=c("vectorclip", "vector", "basic"),
                         A1 = NULL,
                         A1dummy = FALSE,
                         hessian = FALSE,
                         matwt = NULL, new.coef = NULL,
                         saveterms = FALSE,
                         parallel = TRUE
                         ) {
  matrix.action <- match.arg(matrix.action)
  logi.action <- match.arg(logi.action)
  algorithm <- match.arg(algorithm)
  if(reweighting <- !is.null(matwt)) 
    stopifnot(is.numeric(matwt) && is.vector(matwt))
  spill <- spill || spill.vc
  saveterms <- spill && saveterms
  logi <- model$method=="logi"
  asked.parallel <- !missing(parallel)
  
  old.coef <- coef(model)
  use.coef <- if(!is.null(new.coef)) new.coef else old.coef
  p <- length(old.coef)
  if(p == 0) {
    ## this probably can't happen
    if(!spill) return(matrix(, 0, 0)) else return(list())
  }
  pnames <- names(old.coef)
  dnames <- list(pnames, pnames)
  
  internals <- list()
  ##
  sumobj <- summary(model, quick="entries")
  correction <- model$correction
  rbord      <- model$rbord
  R <- reach(model, epsilon=1e-2)
  Q <- quad.ppm(model)
  D <- dummy.ppm(model)
  rho <- model$internal$logistic$rho
  #### If dummy intensity rho is unknown we estimate it
  if(is.null(rho))
     rho <- npoints(D)/(area(D)*markspace.integral(D))
  X <- data.ppm(model)
  Z <- is.data(Q)
  W <- as.owin(model)
  areaW <- if(correction == "border") eroded.areas(W, rbord) else area(W)
  ##
  ## determine which quadrature points contributed to the
  ## sum/integral in the pseudolikelihood
  ## (e.g. some points may be excluded by the border correction)
  okall <- getglmsubset(model)
  ## data only:
  ok <- okall[Z]
  nX <- npoints(X)
  ## conditional intensity lambda(X[i] | X) = lambda(X[i] | X[-i])
  ## data and dummy:
  lamall <- fitted(model, check = FALSE, new.coef = new.coef)
  ## data only:
  lam <- lamall[Z]
  ## sufficient statistic h(X[i] | X) = h(X[i] | X[-i])
  ## data and dummy:
  mall <- model.matrix(model)
  ## save
  if(saveterms) 
    internals <- append(internals,
                        list(mom=mall, lambda=lamall, Z=Z, ok=okall,
                             matwt=matwt))
  if(reweighting) {
    ## each column of the model matrix is multiplied by 'matwt'
    check.nvector(matwt, nrow(mall), things="quadrature points")
    mall.orig <- mall
    mall      <- mall * matwt
  }
  ## subsets of model matrix
  mokall <- mall[okall, , drop=FALSE]
  ## data only:
  m <- mall[Z, , drop=FALSE]
  mok <- m[ok, , drop=FALSE]
  ##
  if(reweighting) {
    ## save unweighted versions
    mokall.orig <- mall.orig[okall, , drop=FALSE]
    m.orig      <- mall.orig[Z, , drop=FALSE]
    mok.orig    <- m.orig[ok, , drop=FALSE]
    ##
    matwtX <- matwt[Z]
  }

  ## ^^^^^^^^^^^^^^^^ First order (sensitivity) matrices A1, S
  
  ## logistic 
  if(logi){
    ## Sensitivity matrix S for logistic case
    Slog <- sumouter(mokall, w = lamall[okall]*rho/(lamall[okall]+rho)^2)
    dimnames(Slog) <- dnames
    ## A1 matrix for logistic case
    A1log <- sumouter(mokall, w = lamall[okall]*rho*rho/(lamall[okall]+rho)^3)
    dimnames(A1log) <- dnames
  }
  ## Sensitivity matrix for MPLE case (= A1)
  if(is.null(A1) || reweighting) {
    if(A1dummy){
      A1 <- sumouter(mokall, w = (lamall * w.quad(Q))[okall])
      if(reweighting)
        gradient <- sumouter(mokall.orig, w=(matwt * lamall * w.quad(Q))[okall])
    } else{
      A1 <- sumouter(mok)
      if(reweighting)
        gradient <- sumouter(mok.orig, w=matwtX)
    }
  } else {
    stopifnot(is.matrix(A1))
    if(!all(dim(A1) == p))
      stop(paste("Matrix A1 has wrong dimensions:",
                 prange(dim(A1)), "!=", prange(c(p, p))))
  }
  dimnames(A1) <- dnames

  ## ^^^^^^^^^^ Second order interaction effects A2, A3

  if(hessian) {
    ## interaction terms suppressed
    A2 <- A3 <- matrix(0, p, p, dimnames=dnames)
    if(logi)
      A2log <- A3log <- matrix(0, p, p, dimnames=dnames)
  } else {
    ## ^^^^^^^^^^^^^^^^^^^^ `parallel' evaluation
    need.loop <- TRUE
    if(parallel) {
      ## compute second order difference
      ##  ddS[i,j,] = h(X[i] | X) - h(X[i] | X[-j])
      ddS <- deltasuffstat(model, restrict=TRUE, force=FALSE)
      if(is.null(ddS)) {
        if(asked.parallel)
          warning("parallel option not available - reverting to loop")
      } else {
        need.loop <- FALSE
        ## rearrange so that
        ##  ddS[ ,i,j] = h(X[i] | X) - h(X[i] | X[-j])
        ddS <- aperm(ddS, c(3,2,1))
        ## now compute sum_{i,j} for i != j
        ## outer(ddS[,i,j], ddS[,j,i])
        ddSok <- ddS[ , ok, ok, drop=FALSE]
        A3 <- sumsymouter(ddSok)
        ## mom.array[ ,i,j] = h(X[i] | X)
        mom.array <- array(t(m), dim=c(p, nX, nX))
        ## momdel[ ,i,j] = h(X[i] | X[-j])
        momdel <- mom.array - ddS
        ## lamdel[i,j] = lambda(X[i] | X[-j])
        lamdel <-
          matrix(lam, nX, nX) * exp(tensor::tensor(-use.coef, ddS, 1, 1))
        ##  pairweight[i,j] = lamdel[i,j]/lambda[i] - 1 
        pairweight <- lamdel / lam - 1
        ## now compute sum_{i,j} for i != j
        ## pairweight[i,j] * outer(momdel[,i,j], momdel[,j,i])
        ## for data points that contributed to the pseudolikelihood
        momdelok <- momdel[ , ok, ok, drop=FALSE]
        A2 <- sumsymouter(momdelok, w=pairweight[ok, ok])
        if(logi){
          ## lam.array[ ,i,j] = lambda(X[i] | X)
          lam.array <- array(lam, c(nX,nX,p))
          lam.array <- aperm(lam.array, c(3,1,2))
          ## lamdel.array[,i,j] = lambda(X[i] | X[-j])
          lamdel.array <- array(lamdel, c(nX,nX,p))
          lamdel.array <- aperm(lamdel.array, c(3,1,2))
          momdellogi <- rho/(lamdel.array+rho)*momdel
          momdellogiok <- momdellogi[ , ok, ok, drop=FALSE]
          A2log <- sumsymouter(momdellogiok, w=pairweight[ok, ok])
          ddSlogi <- rho/(lam.array+rho)*mom.array - momdellogi
          ddSlogiok <- ddSlogi[ , ok, ok, drop=FALSE]
          A3log <- sumsymouter(ddSlogiok)
        }
      }
    }
  
    ## ^^^^^^^^^^^^^^^^^^^^ loop evaluation
    if(need.loop) {
    
      A2 <- A3 <- matrix(0, p, p, dimnames=dnames)
      if(logi)
        A2log <- A3log <- matrix(0, p, p, dimnames=dnames)
    
      if(saveterms) {
        ## *initialise* matrices 
        ##  lamdel[i,j] = lambda(X[i] | X[-j]) = lambda(X[i] | X[-c(i,j)])
        lamdel <- matrix(lam, nX, nX)
        ##  momdel[ ,i,j] = h(X[i] | X[-j]) = h(X[i] | X[-c(i,j)])
        momdel <- array(t(m), dim=c(p, nX, nX))
      }
  
      ## identify close pairs
      if(is.finite(R)) {
        cl <- closepairs(X, R, what="indices")
        I <- cl$i
        J <- cl$j
        if(algorithm == "vectorclip") {
          cl2 <- closepairs(X, 2*R, what="indices")
          I2 <- cl2$i
          J2 <- cl2$j
        }
      } else {
        ## either infinite reach, or something wrong
        IJ <- expand.grid(I=1:nX, J=1:nX)
        IJ <- subset(IJ, I != J)
        I2 <- I <- IJ$I
        J2 <- J <- IJ$J
      }
      ## filter:  I and J must both belong to the nominated subset 
      okIJ <- ok[I] & ok[J]
      I <- I[okIJ]
      J <- J[okIJ]
      ##
      if(length(I) > 0 && length(J) > 0) {
        ## .............. loop over pairs ........................
        ## The following ensures that 'empty' and 'X' have compatible marks 
        empty <- X[integer(0)]
        ## make an empty 'equalpairs' matrix
        nonE <- matrix(, nrow=0, ncol=2)
        ## Run through pairs
        switch(algorithm,
               basic={
                 for(i in unique(I)) {
                   Xi <- X[i]
                   Ji <- unique(J[I==i])
                   if((nJi <- length(Ji)) > 0) {
                     for(k in 1:nJi) {
                       j <- Ji[k]
                       X.ij <- X[-c(i,j)]
                       ## compute conditional intensity
                       ##    lambda(X[j] | X[-i]) = lambda(X[j] | X[-c(i,j)]
                       plamj.i <- predict(model, type="cif",
                                          locations=X[j], X=X.ij,
                                          check = FALSE,
                                          new.coef = new.coef,
                                          sumobj = sumobj, E=nonE)
                       ## corresponding values of sufficient statistic 
                       ##    h(X[j] | X[-i]) = h(X[j] | X[-c(i,j)]
                       pmj.i <- partialModelMatrix(X.ij, X[j], model)[nX-1, ]
                       ## conditional intensity and sufficient statistic
                       ## in reverse order
                       ##    lambda(X[i] | X[-j]) = lambda(X[i] | X[-c(i,j)]
                       plami.j <- predict(model, type="cif",
                                          locations=X[i], X=X.ij,
                                          check = FALSE,
                                          new.coef = new.coef,
                                          sumobj = sumobj, E=nonE)
                       pmi.j <- partialModelMatrix(X.ij, Xi, model)[nX-1, ]
                       ## 
                       if(reweighting) {
                         pmj.i <- pmj.i * matwtX[j]
                         pmi.j <- pmi.j * matwtX[i]
                       }
                       if(saveterms) {
                         lamdel[i,j] <- plami.j
                         momdel[ , i, j] <- pmi.j
                         lamdel[j,i] <- plamj.i
                         momdel[ , j, i] <- pmj.i
                       }
                       ## increment A2, A3
                       wt <- plami.j / lam[i] - 1
                       A2 <- A2 + wt * outer(pmi.j, pmj.i)
                       if(logi)
                         A2log <- A2log +
                           wt * rho/(plami.j+rho) *
                             rho/(plamj.i+rho) * outer(pmi.j, pmj.i)
                       ## delta sufficient statistic
                       ## delta_i h(X[j] | X[-c(i,j)])
                       ## = h(X[j] | X[-j]) - h(X[j] | X[-c(i,j)])
                       ## = h(X[j] | X) - h(X[j] | X[-i])
                       ## delta_j h(X[i] | X[-c(i,j)])
                       ## = h(X[i] | X[-i]) - h(X[i] | X[-c(i,j)])
                       ## = h(X[i] | X) - h(X[i] | X[-j])
                       deltaiSj <- m[j, ] - pmj.i
                       deltajSi <- m[i, ] - pmi.j
                       A3 <- A3 + outer(deltaiSj, deltajSi)
                       if(logi){
                         deltaiSjlog <- rho*(m[j, ]/
                                             (lam[j]+rho) - pmj.i/(plamj.i+rho))
                         deltajSilog <- rho*(m[i, ]/
                                             (lam[i]+rho) - pmi.j/(plami.j+rho))
                         A3log <- A3log + outer(deltaiSjlog, deltajSilog)
                       }
                     }
                   }
                 }
               },
               vector={
                 ## --------- faster algorithm using vector functions --------
                 for(i in unique(I)) {
                   Ji <- unique(J[I==i])
                   nJi <- length(Ji)
                   if(nJi > 0) {
                     Xi <- X[i]
                     ## neighbours of X[i]
                     XJi <- X[Ji]
                     ## all points other than X[i]
                     X.i <- X[-i]
                     ## index of XJi in X.i
                     J.i <- Ji - (Ji > i)
                     ## equalpairs matrix
                     E.i <- cbind(J.i, seq_len(nJi))
                     ## compute conditional intensity
                     ##   lambda(X[j] | X[-i]) = lambda(X[j] | X[-c(i,j)]
                     ## for all j
                     plamj <- predict(model, type="cif",
                                      locations=XJi, X=X.i,
                                      check = FALSE,
                                      new.coef = new.coef,
                                      sumobj=sumobj, E=E.i)
                     ## corresponding values of sufficient statistic 
                     ##    h(X[j] | X[-i]) = h(X[j] | X[-c(i,j)]
                     ## for all j
                     pmj <-
                       partialModelMatrix(X.i, empty, model)[J.i, , drop=FALSE]
                     ##
                     ## conditional intensity & sufficient statistic
                     ## in reverse order
                     ##    lambda(X[i] | X[-j]) = lambda(X[i] | X[-c(i,j)]
                     ## for all j
                     plami <- numeric(nJi)
                     pmi <- matrix(, nJi, p)
                     for(k in 1:nJi) {
                       j <- Ji[k]
                       X.ij <- X[-c(i,j)]
                       plami[k] <- predict(model, type="cif",
                                           locations=Xi, X=X.ij,
                                           check = FALSE,
                                           new.coef = new.coef,
                                           sumobj = sumobj, E=nonE)
                       pmi[k, ] <- partialModelMatrix(X.ij, Xi, model)[nX-1, ]
                     }
                     ##
                     if(reweighting) {
                       pmj <- pmj * matwtX[Ji]
                       pmi <- pmi * matwtX[i]
                     }
                     if(saveterms) {
                       lamdel[Ji, i] <- plamj
                       momdel[ , Ji, i] <- t(pmj)
                       lamdel[i,Ji] <- plami
                       momdel[ , i, Ji] <- t(pmi)
                     }
                     ## increment A2, A3
                     wt <- plami / lam[i] - 1
                     for(k in 1:nJi) {
                       j <- Ji[k]
                       A2 <- A2 + wt[k] * outer(pmi[k,], pmj[k,])
                       if(logi)
                         A2log <- A2log + wt[k] * rho/(plami[k]+rho) *
                           rho/(plamj[k]+rho) * outer(pmi[k,], pmj[k,])
                       ## delta sufficient statistic
                       ## delta_i h(X[j] | X[-c(i,j)])
                       ## = h(X[j] | X[-j]) - h(X[j] | X[-c(i,j)])
                       ## = h(X[j] | X) - h(X[j] | X[-i])
                       ## delta_j h(X[i] | X[-c(i,j)])
                       ## = h(X[i] | X[-i]) - h(X[i] | X[-c(i,j)])
                       ## = h(X[i] | X) - h(X[i] | X[-j])
                       deltaiSj <- m[j, ] - pmj[k,]
                       deltajSi <- m[i, ] - pmi[k,]
                       A3 <- A3 + outer(deltaiSj, deltajSi)
                       if(logi){
                         deltaiSjlog <- rho*(m[j, ]/(lam[j]+rho) -
                                             pmj[k,]/(plamj[k]+rho))
                         deltajSilog <- rho*(m[i, ]/(lam[i]+rho) -
                                             pmi[k,]/(plami[k]+rho))
                         A3log <- A3log + outer(deltaiSjlog, deltajSilog)
                       }
                     }
                   }
                 }
               },
               vectorclip={
                 ## --------- faster version of 'vector' algorithm
                 ## --------  by removing non-interacting points of X
                 for(i in unique(I)) {
                   ## all points within 2R
                   J2i <- unique(J2[I2==i])
                   ## all points within R
                   Ji  <- unique(J[I==i])
                   nJi <- length(Ji)
                   if(nJi > 0) {
                     Xi <- X[i]
                     ## neighbours of X[i]
                     XJi <- X[Ji]
                     ## replace X[-i] by X[-i] \cap b(0, 2R)
                     X.i <- X[J2i]
                     nX.i <- length(J2i)
                     ## index of XJi in X.i
                     J.i <- match(Ji, J2i)
                     if(any(is.na(J.i)))
                       stop("Internal error: Ji not a subset of J2i")
                     ## equalpairs matrix
                     E.i <- cbind(J.i, seq_len(nJi))
                     ## compute conditional intensity
                     ##   lambda(X[j] | X[-i]) = lambda(X[j] | X[-c(i,j)]
                     ## for all j
                     plamj <- predict(model, type="cif",
                                      locations=XJi, X=X.i,
                                      check = FALSE,
                                      new.coef = new.coef,
                                      sumobj = sumobj, E=E.i)
                     ## corresponding values of sufficient statistic 
                     ##    h(X[j] | X[-i]) = h(X[j] | X[-c(i,j)]
                     ## for all j
                     pmj <-
                       partialModelMatrix(X.i, empty, model)[J.i, , drop=FALSE]
                     ##
                     ## conditional intensity & sufficient statistic
                     ##  in reverse order
                     ##    lambda(X[i] | X[-j]) = lambda(X[i] | X[-c(i,j)]
                     ## for all j
                     plami <- numeric(nJi)
                     pmi <- matrix(, nJi, p)
                     for(k in 1:nJi) {
                       j <- Ji[k]
                       ## X.ij <- X[-c(i,j)]
                       X.ij <- X.i[-J.i[k]]
                       plami[k] <- predict(model, type="cif",
                                           locations=Xi, X=X.ij,
                                           check = FALSE,
                                           new.coef = new.coef,
                                           sumobj = sumobj, E=nonE)
                       pmi[k, ] <- partialModelMatrix(X.ij, Xi, model)[nX.i, ]
                     }
                     ##
                     if(reweighting) {
                       pmj <- pmj * matwtX[Ji]
                       pmi <- pmi * matwtX[i]
                     }
                     if(saveterms) {
                       lamdel[Ji, i] <- plamj
                       momdel[ , Ji, i] <- t(pmj)
                       lamdel[i,Ji] <- plami
                       momdel[ , i, Ji] <- t(pmi)
                     }
                     ## increment A2, A3
                     wt <- plami / lam[i] - 1
                     for(k in 1:nJi) {
                       j <- Ji[k]
                       A2 <- A2 + wt[k] * outer(pmi[k,], pmj[k,])
                       if(logi)
                         A2log <- A2log + wt[k] * rho/(plami[k]+rho) *
                           rho/(plamj[k]+rho) * outer(pmi[k,], pmj[k,])
                       ## delta sufficient statistic
                       ## delta_i h(X[j] | X[-c(i,j)])
                       ## = h(X[j] | X[-j]) - h(X[j] | X[-c(i,j)])
                       ## = h(X[j] | X) - h(X[j] | X[-i])
                       ## delta_j h(X[i] | X[-c(i,j)])
                       ## = h(X[i] | X[-i]) - h(X[i] | X[-c(i,j)])
                       ## = h(X[i] | X) - h(X[i] | X[-j])
                       deltaiSj <- m[j, ] - pmj[k,]
                       deltajSi <- m[i, ] - pmi[k,]
                       A3 <- A3 + outer(deltaiSj, deltajSi)
                       if(logi){
                         deltaiSjlog <- rho*(m[j, ]/(lam[j]+rho) -
                                             pmj[k,]/(plamj[k]+rho))
                         deltajSilog <- rho*(m[i, ]/(lam[i]+rho) -
                                             pmi[k,]/(plami[k]+rho))
                         A3log <- A3log + outer(deltaiSjlog, deltajSilog)
                       }
                     }
                   }
                 }
               })
      }
    }
    ## ......... end of loop computation ...............
  }

  #### Matrix Sigma 
  Sigma <- A1+A2+A3
  
  if(spill) {
    ## save internal data (with matrices unnormalised) 
    internals <-
      c(internals,
        list(A1=A1, A2=A2, A3=A3, Sigma=Sigma, areaW=areaW),
        if(logi)
           list(A1log=A1log, A2log=A2log, A3log=A3log, Slog=Slog) else NULL,
        if(reweighting) list(gradient=gradient) else NULL,
        if(saveterms) list(lamdel=lamdel, momdel=momdel) else NULL)
    ## return internal data if no further calculation needed
    if(!spill.vc && !logi)
      return(internals)
  }
    
  ## ........... calculate variance/covariance matrix for MPL .........

  if(!reweighting) {
    ## Normalise
    A1 <- A1/areaW
    Sigma <- Sigma/areaW
    ## Enforce exact symmetry 
    A1 <- (A1 + t(A1))/2
    Sigma <- (Sigma + t(Sigma))/2
    ## calculate inverse negative Hessian
    U <- checksolve(A1, matrix.action, , "variance")
  } else {
    ## Normalise
    gradient <- gradient/areaW
    Sigma <- Sigma/areaW
    ## Enforce exact symmetry
    gradient <- (gradient + t(gradient))/2
    Sigma <- (Sigma + t(Sigma))/2
    ## calculate inverse negative Hessian
    U <- checksolve(gradient, matrix.action, , "variance")
  }
  
  ## compute variance-covariance
  vc.mpl <- if(is.null(U)) matrix(NA, p, p) else 
              U %*% Sigma %*% U / areaW
  dimnames(vc.mpl) <- dnames

  ## return variance-covariance matrix, if model was fitted by MPL
  if(!logi) {
    if(spill.vc) return(list(varcov=vc.mpl, internals=internals))
    return(vc.mpl)
  }
  
  ###### Everything below is only computed for logistic fits #######

  ## Matrix Sigma1log (A1log+A2log+A3log):
  Sigma1log <- A1log+A2log+A3log
  ## Resolving the dummy process type
  how <- model$internal$logistic$how
  if(how %in% c("given", "grid", "transgrid")){
    whinge <- paste("vcov is not implemented for dummy type", sQuote(how))
    if(logi.action=="fatal")
      stop(whinge)
    how <- if(how=="given") "poisson" else "stratrand"
    if(logi.action=="warn")
      warning(paste(whinge,"- using", sQuote(how), "formula"), call.=FALSE)
  }
  ## Matrix Sigma2log (depends on dummy process type)
  switch(how,
         poisson={
           Sigma2log <- sumouter(mokall, w = lamall[okall]*lamall[okall]*rho/(lamall[okall]+rho)^3)
         },
         binomial={
           Sigma2log <- sumouter(mokall, w = lamall[okall]*lamall[okall]*rho/(lamall[okall]+rho)^3)
           A1vec <- t(mokall) %*% (rho*lamall[okall]/(lamall[okall]+rho)^2)
           Sigma2log <- Sigma2log - A1vec%*%t(A1vec)/rho*1/sum(1/(lamall[okall]+rho))
         },
         stratrand={
           ## Dirty way of refitting model with new dummy pattern (should probably be done using call, eval, envir, etc.):
           ## Changed by ER 2013/06/14 to use the new quadscheme.logi
           ## D2 <- logi.dummy(X = X, type = "stratrand", nd = model$internal$logistic$args)
           ## Q2 <- quad(data=X, dummy=D2)
           ## Q2$dummy$Dinfo <- D2$Dinfo
           Q2 <- quadscheme.logi(data=X, dummytype = "stratrand",
                                 nd = model$internal$logistic$nd)
           D2 <- Q2$dummy
           Q2$dummy$Dinfo <- D2$Dinfo
           Z2 <- is.data(Q2)
           arglist <- list(Q=Q2, trend=model$trend, interaction = model$interaction, method = model$method,
                           correction = model$correction, rbord = model$rbord, covariates = model$covariates)
           arglist <- append(arglist, model$internal$logistic$extraargs)
           model2 <- do.call(ppm, args = arglist)

           ## New cif
           lamall2 <- fitted(model2, check = FALSE, new.coef = new.coef)
           ## New model matrix
           mall2 <- model.matrix(model2)
           okall2 <- getglmsubset(model2)

           ## index vectors of stratrand cell indices of dummy points 
           inD <- model$internal$logistic$inD
           inD2 <- model2$internal$logistic$inD

           ## Dummy points inside eroded window (for border correction)
           if(is.finite(R) && (correction == "border")){
             ii <- (bdist.points(D) >= R)
             ii2 <- (bdist.points(D2) >= R)
           } else{
             ii <- rep.int(TRUE, npoints(D))
             ii2 <- rep.int(TRUE, npoints(D2))
           }
           ## OK points of dummy pattern 1 with a valid point of dummy pattern 2 in same stratrand cell (and vice versa)
           okdum <- okall[!Z]
           okdum2 <- okall2[!Z2]
           ok1 <- okdum & ii & is.element(inD, inD2[okdum2 & ii2])
           ok2 <- okdum2 & ii2 & is.element(inD2, inD[okdum & ii])
           ## ok1 <- okdum & okdum2 & ii & is.element(inD, inD2[ii2])
           ## ok2 <- okdum2 & okdum1 & ii2 & is.element(inD2, inD[ii])
           ## ok1 <- ii & is.element(inD, inD2[ii2])
           ## ok2 <- ii2 & is.element(inD2, inD[ii])

           ## cif and suff. stat. for valid points in dummy patterns 1 and 2
           lamdum <- lamall[!Z][ok1]
           lamdum2 <- lamall2[!Z2][ok2]
           mdum <- mall[!Z,,drop=FALSE][ok1,]
           mdum2 <- mall2[!Z2,,drop=FALSE][ok2,]

           ## finally calculation of Sigma2
           wlam <- mdum * rho*lamdum/(lamdum+rho)
           wlam2 <- mdum2 * rho*lamdum2/(lamdum2+rho)
           Sigma2log <- t(wlam-wlam2)%*%(wlam-wlam2)/(2*rho*rho)
         },
         stop("sorry - unrecognized dummy process in logistic fit")
         )
  ## Attaching to Sigma2log calculated above
  dimnames(Sigma2log) <- dnames

  
  if(spill) {
    ## return internal data only (with matrices unnormalised)
    internals <- c(internals, 
                   list(Sigma1log=Sigma1log, Sigma2log=Sigma2log, mple=vc.mpl))
    if(!spill.vc)
      return(internals)
  }

  ## .. Calculate variance-covariance matrix for logistic fit ...........
  ## normalise
  Slog <- Slog/areaW
  Sigma1log <- Sigma1log/areaW
  Sigma2log <- Sigma2log/areaW
  ## evaluate
  Ulog <- checksolve(Slog, matrix.action, , "variance")
  vc.logi <- if(is.null(Ulog)) matrix(NA, p, p) else 
             Ulog %*% (Sigma1log+Sigma2log) %*% Ulog / areaW
  dimnames(vc.logi) <- dnames
  ##
  if(spill.vc) return(list(varcov=vc.logi, internals=internals))
  return(vc.logi)
}

## vcalcGibbs from Ege Rubak and J-F Coeurjolly
## 2013/06/14, modified by Ege to handle logistic case as well

vcalcGibbsSpecial <- function(fit, ...,
                              spill=FALSE,
                              spill.vc=FALSE,
                              special.alg = TRUE,
                              matrix.action=c("warn", "fatal", "silent"),
                              logi.action=c("warn", "fatal", "silent")) {
  matrix.action <- match.arg(matrix.action)
  logi.action <- match.arg(logi.action)
  spill <- spill || spill.vc
  
  ## Interaction name:
  iname <- fit$interaction$name
  
  ## Does the model have marks which are in the trend?
  marx <- is.marked(fit) && ("marks" %in% variablesinformula(fit$trend))

  ## The full data and window:
  Xplus <- data.ppm(fit)
  Wplus <- as.owin(Xplus)

  ## Fitted parameters and the parameter dimension p (later consiting of p1 trend param. and p2 interaction param.):
  theta <- coef(fit)
  p <- length(theta)

  ## Number of points:
  n <- npoints(Xplus)

  ## Using the faster algorithms for special cases
  if(special.alg && fit$method != "logi"){
    param <- coef(fit)
    switch(iname,
      "Strauss process"={
        ## Only implemented for non-marked case:
        if(!marx)
	  return(vcovPairPiece(Xplus,
                               reach(fit$interaction),
                               exp(coef(fit)[2]),
                               matrix.action,
                               spill=spill,
                               spill.vc=spill.vc))
      },
           
      "Piecewise constant pairwise interaction process"={
        ## Only implemented for non-marked case:
        if(!marx)
          return(vcovPairPiece(Xplus,
                               fit$interaction$par$r,
                               exp(coef(fit)[-1]),
                               matrix.action,
                               spill=spill,
                               spill.vc=spill.vc))
      },

      "Multitype Strauss process"={
	matR <- fit$interaction$par$radii
        R <- c(matR[1,1], matR[1,2], matR[2,2])
        ## Only implemented for 2 types with equal interaction range:
        if(ncol(matR)==2 && marx){
          n <- length(theta)
          res <- vcovMultiStrauss(Xplus, R, exp(theta[c(n-2,n-1,n)]),
                                  matrix.action,spill=spill,spill.vc=spill.vc)
          if(!spill) {
            res <- contrastmatrix(res, 2)
            dimnames(res) <- list(names(theta), names(theta))
          }
          return(res)
        }
      }
    )
  }
  
  ## Matrix specifying equal points in the two patterns in the call to eval below:
  E <- matrix(rep.int(1:n, 2), ncol = 2)

  ## Eval. the interaction potential difference at all points (internal spatstat function):
#  V1 <- fit$interaction$family$eval(Xplus, Xplus, E, fit$interaction$pot, fit$interaction$par, fit$correction)
  oldopt <- NULL
  if(fit$interaction$family$name=="pairwise"){
      oldopt <- spatstat.options(fasteval = "off")
  }
  V1 <- evalInteraction(Xplus, Xplus, E, as.interact(fit), fit$correction)
  spatstat.options(oldopt)

  ## Calculate parameter dimensions and correct the contrast type parameters:
  p2 <- ncol(V1)
  p1 <- p-p2
  if(p1>1)
    theta[2:p1] <- theta[2:p1] + theta[1]
  ## V1 <- evalInteraction(Q, Xplus, union.quad(Q), fit$interaction, fit$correction)
  POT <- attr(V1, "POT")
  attr(V1, "POT") <- NULL
  ## Adding the constant potential as first column (one column per type for multitype):
  if(!marx){
    V1 <- cbind(1, V1)
    colnames(V1) <- names(theta)
  }
  else{
    lev <- levels(marks(Xplus))
    ## Indicator matrix for mark type attached to V1:
    tmp <- matrix(marks(Xplus), nrow(V1), p1)==matrix(lev, nrow(V1), p-ncol(V1), byrow=TRUE)
    colnames(tmp) <- lev
    V1 <- cbind(tmp,V1)
  }

  ## Matrices for differences of potentials:
  E <- matrix(rep.int(1:(n-1), 2), ncol = 2)
  dV <- V2 <- array(0,dim=c(n,n,p))

  for(k in 1:p1){
    V2[,,k] <- matrix(V1[,k], n, n, byrow = FALSE)
  }
  for(k in (p1+1):p){
    diag(V2[,,k]) <- V1[,k]
  }
  for(j in 1:n){
    ## Fast evaluation for pairwise interaction processes:
    if(fit$interaction$family$name=="pairwise" && !is.null(POT)){
      V2[-j,j,-(1:p1)] <- V1[-j,-(1:p1)]-POT[-j,j,]
    }
    else{
      V2[-j,j,-(1:p1)] <- fit$interaction$family$eval(Xplus[-j], Xplus[-j], E, fit$interaction$pot, fit$interaction$par, fit$correction)
      ## Q <- quadscheme(Xplus[-j],emptyppp)
      ## V2[-j,j,-1] <- evalInteraction(Q, Xplus[-j], Xplus[-j], fit$interaction, fit$correction)
    }
    for(k in 1:p){
      dV[,j,k] <- V1[,k] - V2[,j,k]
    }
  }
  ## Ratio of first and second order Papangelou - 1:
  frac <- 0*dV[,,1]
  for(k in (p1+1):p){
    frac <- frac + dV[,,k]*theta[k]
  }
  frac <- exp(-frac)-1

  ## In the rest we restrict attention to points in the interior:
  
  ## The interaction range:
  R <- reach(fit$interaction)

  ## The reduced window, area and point pattern:
  W<-erosion.owin(Wplus,R)
  areaW <- area(W)

  ## Interior points determined by bdist.points:
  IntPoints <- bdist.points(Xplus)>=R  
  X <- Xplus[IntPoints]
  
  ## Making a logical matrix, I, indicating R-close pairs which are in the interior:
  D <- pairdist(Xplus)
  diag(D) <- Inf
  I <- (D<=R) & outer(IntPoints,IntPoints, "&")
  
  ## Matrix A1:
  A1 <- t(V1[IntPoints,])%*%V1[IntPoints,]

  ## Matrix A2:
  A2 <- matrix(0,p,p)
  for(k in 1:p){
    for(l in k:p){
      A2[k,l] <- A2[l,k] <- sum(I*V2[,,k]*frac*t(V2[,,l]))
    }
  }
  
  ## Matrix A3:
  A3 <- matrix(0,p,p)
  for(k in 1:p){
    for(l in k:p){
      A3[k,l] <- A3[l,k] <- sum(I*dV[,,k]*t(dV[,,l]))
    }
  }

  ## Matrix Sigma (A1+A2+A3):
  Sigma<-A1+A2+A3

  if(spill) {
    # save internal data (with matrices unnormalised)
    dimnames(A1) <- dimnames(A2) <-
      dimnames(A3) <- list(names(theta), names(theta))
    internals <- list(A1=A1, A2=A2, A3=A3, Sigma=Sigma, areaW=areaW)
    # return internal data, if model fitted by MPL
    if(!spill.vc && fit$method != "logi")
      return(internals)
  }

  # ......... Calculate variance-covariance matrix for MPL ........
  
  # normalise
  A1 <- A1/areaW
  Sigma <- Sigma/areaW
  # evaluate
  U <- checksolve(A1, matrix.action, , "variance")
  vc.mpl <- if(is.null(U)) matrix(NA, p, p) else U %*% Sigma %*% U / areaW
  ## Convert to treatment contrasts
  if(marx)
    vc.mpl <- contrastmatrix(vc.mpl, p1)
  dimnames(vc.mpl) <- list(names(theta), names(theta))
  
  # Return result for standard ppm method:
  if(fit$method!="logi") {
    if(spill.vc) return(list(varcov=vc.mpl, internals=internals))
    return(vc.mpl)
  }
  
  ########################################################################
  ###### The remainder is only executed when the method is logistic ######
  ########################################################################

  ### Most of this is copy/pasted from vcalcGibbsGeneral
  correction <- fit$correction
  Q <- quad.ppm(fit)
  D <- dummy.ppm(fit)
  rho <- fit$internal$logistic$rho
  ## If dummy intensity rho is unknown we estimate it
  if(is.null(rho))
     rho <- npoints(D)/(area(D)*markspace.integral(D))
  X <- data.ppm(fit)
  Z <- is.data(Q)

  # determine which data points entered into the sum in the pseudolikelihood
  # (border correction, nonzero cif)
  # data and dummy:
  okall <- getglmsubset(fit)
  ## # data only:
  ## ok <- okall[Z]

  # conditional intensity lambda(X[i] | X) = lambda(X[i] | X[-i])
  # data and dummy:
  lamall <- fitted(fit, check = FALSE)
  ## # data only:
  ## lam <- lamall[Z]

  # sufficient statistic h(X[i] | X) = h(X[i] | X[-i])
  # data and dummy:
  mall <- model.matrix(fit)
  mokall <- mall[okall, , drop=FALSE]
  ## # data only:
  ## m <- mall[Z, , drop=FALSE]
  ## mok <- m[ok, , drop=FALSE]

  # Sensitivity matrix S and A1 matrix for logistic case
  Slog <- sumouter(mokall, w = lamall[okall]*rho/(lamall[okall]+rho)^2)
  A1log <- sumouter(mokall, w = lamall[okall]*rho*rho/(lamall[okall]+rho)^3)

  ## Define W1, W2 and dW for the logistic method based on V1, V2 and dV (frac is unchanged)
  lambda1 <- exp(.rowSums(matrix(theta,n,p,byrow=TRUE)*V1, n, p))
  W1 <- V1*rho/(lambda1+rho)
  lambda2 <- exp(apply(array(rep(theta,each=n*n),dim=c(n,n,p))*V2, c(1,2), sum))
  W2 <- V2
  dW <- dV
  for(k in 1:p){
    W2[,,k] <- V2[,,k] * rho/(lambda2+rho)
    for(j in 1:n){
      dW[,j,k] <- W1[,k] - W2[,j,k]
    }
  }
  ## Matrices A2log and A3log for the first component Sigma1log of the variance:
  A2log <- A3log <- matrix(0,p,p)
  for(k in 1:p){
    for(l in k:p){
      A2log[k,l] <- A2log[l,k] <- sum(I*W2[,,k]*frac*t(W2[,,l]))
      A3log[k,l] <- A3log[l,k] <- sum(I*dW[,,k]*t(dW[,,l]))
    }
  }
  A2log <- A2log
  A3log <- A3log
  
  ## First variance component Sigma1log (A1log+A2log+A3log):
  Sigma1log <- A1log+A2log+A3log

  ## Resolving the dummy process type
  how <- fit$internal$logistic$how
  if(how %in% c("given", "grid", "transgrid")){
    whinge <- paste("vcov is not implemented for dummy type", sQuote(how))
    if(logi.action=="fatal")
      stop(whinge)
    how <- if(how=="given") "poisson" else "stratrand"
    if(logi.action=="warn")
      warning(paste(whinge,"- using", sQuote(how), "formula"), call.=FALSE)
  }

  ## Matrix Sigma2log (depends on dummy process type)
  switch(how,
         poisson={
           Sigma2log <- sumouter(mokall, w = lamall[okall]*lamall[okall]*rho/(lamall[okall]+rho)^3)
         },
         binomial={
           Sigma2log <- sumouter(mokall, w = lamall[okall]*lamall[okall]*rho/(lamall[okall]+rho)^3)
           A1vec <- t(mokall) %*% (rho*lamall[okall]/(lamall[okall]+rho)^2)
           Sigma2log <- Sigma2log - A1vec%*%t(A1vec)/rho*1/sum(1/(lamall[okall]+rho))
         },
         stratrand={
           ### Dirty way of refitting model with new dummy pattern (should probably be done using call, eval, envir, etc.):
           ## D2 <- logi.dummy(X = X, type = "stratrand", nd = model$internal$logistic$args)
           ## Q2 <- quad(data=X, dummy=D2)
           ## Q2$dummy$Dinfo <- D2$Dinfo
           Q2 <- quadscheme.logi(data=X, dummytype = "stratrand", nd = fit$internal$logistic$nd)
           D2 <- Q2$dummy
           Z2 <- is.data(Q2)
           arglist <- list(Q=Q2, trend=fit$trend, interaction = fit$interaction, method = fit$method,
                           correction = fit$correction, rbord = fit$rbord, covariates = fit$covariates)
           arglist <- append(arglist, fit$internal$logistic$extraargs)
           fit2 <- do.call(ppm, args = arglist)

           ## New cif
           lamall2 <- fitted(fit2, check=FALSE)
           ## New model matrix
           mall2 <- model.matrix(fit2)
           okall2 <- getglmsubset(fit2)

           # index vectors of stratrand cell indices of dummy points 
           inD <- fit$internal$logistic$inD
           inD2 <- fit2$internal$logistic$inD

           # Dummy points inside eroded window (for border correction)
           if(is.finite(R) && (correction == "border")){
             ii <- inside.owin(D, w = W)
             ii2 <- inside.owin(D2, w = W)
           } else{
             ii <- rep.int(TRUE, npoints(D))
             ii2 <- rep.int(TRUE, npoints(D2))
           }
           # OK points of dummy pattern 1 with a valid point of dummy pattern 2 in same stratrand cell (and vice versa)
           okdum <- okall[!Z]
           okdum2 <- okall2[!Z2]
           ok1 <- okdum & ii & is.element(inD, inD2[okdum2 & ii2])
           ok2 <- okdum2 & ii2 & is.element(inD2, inD[okdum & ii])
           ## ok1 <- okdum & okdum2 & ii & is.element(inD, inD2[ii2])
           ## ok2 <- okdum2 & okdum1 & ii2 & is.element(inD2, inD[ii])
           ## ok1 <- ii & is.element(inD, inD2[ii2])
           ## ok2 <- ii2 & is.element(inD2, inD[ii])

           # cif and suff. stat. for valid points in dummy patterns 1 and 2
           lamdum <- lamall[!Z][ok1]
           lamdum2 <- lamall2[!Z2][ok2]
           mdum <- mall[!Z,][ok1,]
           mdum2 <- mall2[!Z2,][ok2,]

           # finally calculation of Sigma2
           wlam <- mdum * rho*lamdum/(lamdum+rho)
           wlam2 <- mdum2 * rho*lamdum2/(lamdum2+rho)
           Sigma2log <- t(wlam-wlam2)%*%(wlam-wlam2)/(2*rho*rho)
         },
         stop("sorry - unrecognized dummy process in logistic fit")
         )


  if(spill) {
    ## Attach dimnames to all matrices
    dimnames(Sigma2log) <- dimnames(Slog) <-
      dimnames(Sigma1log) <- dimnames(A1log) <-
        dimnames(A2log) <- dimnames(A3log) <-
          list(names(theta),names(theta))
    # return internal data (with matrices unnormalised)
    internals <- c(internals,
                   list(A1log=A1log, A2log=A2log, A3log=A3log, Slog=Slog,
                        Sigma1log=Sigma1log, Sigma2log=Sigma2log, mple=vc.mpl))
    if(!spill.vc)
      return(internals)
  }

  # ....... Compute variance-covariance for logistic fit .............
  # Normalise
  Slog <- Slog/areaW
  Sigma1log <- Sigma1log/areaW
  Sigma2log <- Sigma2log/areaW
  ## Finally the result is calculated:
  Ulog <- checksolve(Slog, matrix.action, , "variance")
  vc.logi <- if(is.null(Ulog)) matrix(NA, p, p) else 
             Ulog %*% (Sigma1log+Sigma2log) %*% Ulog / areaW
  #
  dimnames(vc.logi) <- list(names(theta), names(theta))
  if(spill.vc) return(list(varcov=vc.logi, internals=internals))
  return(vc.logi)
}

vcovPairPiece <- function(Xplus, R, Gam, matrix.action,
                          spill=FALSE, spill.vc=FALSE){
  ## R is  the  vector of breaks (R[length(R)]= range of the pp.
  ## Gam is the vector of weights
  Rmax <- R[length(R)]
  
  ## Xplus : point process observed in W+R
  ## Extracting the window and calculating area:
  Wplus<-as.owin(Xplus)
  W<-erosion.owin(Wplus,Rmax)
  areaW <- area(W)

  ## Interior points determined by bdist.points:
  IntPoints <- bdist.points(Xplus)>=Rmax
  X <- Xplus[IntPoints]

  nX <- npoints(X)
  nXplus <- npoints(Xplus)
  ## Matrix D with pairwise distances between points and infinite distance
  ## between a point and itself:
  
  Dplus<-pairdist(Xplus)
  D <- pairdist(X)
  diag(D) <- diag(Dplus) <- Inf
  ## logical matrix, I, indicating R-close pairs:
  p<-length(R)
  Tplus<-T<-matrix(0,X$n,p)
  I<-Iplus<-list()
  for (i in 1:p){
     if (i==1){
	Iplus[[1]]<- Dplus <=R[1]
	I[[1]] <- D<=R[1]
     } else {
	Iplus[[i]]<- ((Dplus>R[i-1]) & (Dplus <=R[i]))
	I[[i]] <- ((D>R[i-1]) & (D <=R[i]))
     }
     ## Vector T with the number of $R$-close neighbours to each point:
     Tplus[,i]<-  .colSums(Iplus[[i]], nXplus, nXplus)[IntPoints]
     T[,i] <-  .colSums(I[[i]], nX, nX)
  }
  ## Matrices A1, A2 and A3 are initialized to zero:
  A1 <- A2 <- A3 <- matrix(0,p+1,p+1)
  ## A1 and A3:
  A1[1,1] <- npoints(X)
  
  for (j in (2:(p+1))){
    A1[1,j]<-A1[j,1]<-sum(Tplus[,j-1])
    A3[j,j]<-sum(T[,j-1])
    for (k in (2:(p+1))){
      A1[j,k]<-sum(Tplus[,j-1] * Tplus[,k-1])
    }
  }
  ## A2:
  for (j in (2:(p+1))){
    A2[1,1]<-A2[1,1]+(Gam[j-1]^(-1)-1)*sum(T[,j-1])
    for (l in (2:(p+1))){
      if (l==j) vj<-Tplus[,j-1]-1 else vj<-Tplus[,j-1]
	A2[1,j]<-A2[1,j]+(Gam[l-1]^(-1)-1)*sum(T[,l-1]*(vj) )
    }
    A2[j,1]<-A2[1,j]
    for (k in (2:(p+1))){
      for (l in (2:(p+1))){
	if (l==j) vj<-Tplus[,j-1]-1 else vj<-Tplus[,j-1]
	if (l==k) vk<-Tplus[,k-1]-1 else vk<-Tplus[,k-1]

	A2[j,k]<-A2[j,k]+ (Gam[l-1]^(-1)-1)*sum(I[[l-1]]*outer(vj,vk))
      }
    }

  }

  Sigma<-A1+A2+A3

  nam <- c("(Intercept)", names(Gam))
  dnam <- list(nam, nam)
  
  if(spill) {
    # return internal data (with matrices unnormalised)
    dimnames(A1) <- dimnames(A2) <- dimnames(A3) <- dimnames(Sigma) <- dnam
    internals <- list(A1=A1, A2=A2, A3=A3, Sigma=Sigma)
    if(!spill.vc) return(internals)
  }
           
  ## Calculate variance-covariance
  # Normalise:
  A1    <- A1/areaW
  Sigma <- Sigma/areaW
  U <- checksolve(A1, matrix.action, , "variance")
  mat <- if(is.null(U)) matrix(NA, length(nam), length(nam)) else U%*%Sigma%*%U / areaW
  dimnames(mat) <- dnam

  if(spill.vc) return(list(varcov=mat, internals=internals))
  return(mat)
}

vcovMultiStrauss <- function(Xplus, vecR, vecg, matrix.action,
                             spill=FALSE, spill.vc=FALSE){
  ## Xplus : marked Strauss point process 
  ## with two types 
  ## observed in W+R (R=max(R11,R12,R22))

  ## vecg =  estimated parameters of interaction parameters
  ##	    ordered as the output of ppm, i.e. vecg=(g11,g12,g22)	
  ## vecR = range for the diff. strauss ordered a vecg(R11,R12,R22)

  R <- max(vecR)
  R11<-vecR[1];R12<-vecR[2];R22<-vecR[3]
  ## Extracting the window and calculating area:
  Wplus<-as.owin(Xplus)
  W<-erosion.owin(Wplus,R)
  areaW <- area(W)
  X1plus<-Xplus[Xplus$marks==levels(Xplus$marks)[1]]
  X2plus<-Xplus[Xplus$marks==levels(Xplus$marks)[2]]

  ## Interior points determined by bdist.points:
  IntPoints1 <- bdist.points(X1plus)>=R
  IntPoints2 <- bdist.points(X2plus)>=R
  X1 <- X1plus[IntPoints1]
  X2 <- X2plus[IntPoints2]

  nX1 <- npoints(X1)
  nX2 <- npoints(X2)
  nX1plus <- npoints(X1plus)
  nX2plus <- npoints(X2plus)
  
  ## Matrix D with pairwise distances between points and infinite distance
  ## between a point and itself:

  D1plus<-pairdist(X1plus)
  D1 <- pairdist(X1)
  diag(D1) <- diag(D1plus) <- Inf
  
  D2plus<-pairdist(X2plus)
  D2 <- pairdist(X2)
  diag(D2) <- diag(D2plus) <- Inf
  
  D12plus<-crossdist(X1,X2plus)  
  T12plus<-  .rowSums(D12plus<=R12, nX1, nX2plus)
  D21plus<-crossdist(X2,X1plus) 
  T21plus<-  .rowSums(D21plus<=R12, nX2, nX1plus)
  
  I12<-crossdist(X1,X2)<=R12
  I21<-crossdist(X2,X1)<=R12
  T12<-   .rowSums(I12, nX1, nX2)
  T21<-   .rowSums(I21, nX2, nX1)
  ## logical matrix, I, indicating R-close pairs:
  I1plus<- D1plus <=R11
  I1 <- D1<=R11
  I2plus<- D2plus <=R22
  I2 <- D2<=R22
  ## Vector T with the number of $R$-close neighbours to each point:
  T1plus<-  .colSums(I1plus, nX1plus, nX1plus)[IntPoints1]
  T1 <-     .colSums(I1,     nX1,     nX1)
  T2plus<-  .colSums(I2plus, nX2plus, nX2plus)[IntPoints2]
  T2 <-     .colSums(I2,     nX2,     nX2)

  ## Matrices A1, A2 and A3 are initialized to zero:
  A1 <- A2 <- A3 <- matrix(0,5,5)
  ## A1 is filled:
  A1[1,1]<-npoints(X1)
  A1[1,3]<-A1[3,1]<-sum(T1plus)
  A1[1,4]<-A1[4,1]<-sum(T12plus)
  A1[2,2]<-npoints(X2)
  A1[2,5]<-A1[5,2]<-sum(T2plus)
  A1[2,4]<-A1[4,2]<-sum(T21plus)
  A1[3,3]<-sum(T1plus*T1plus)
  A1[3,4]<-A1[4,3]<-sum(T1plus*T12plus)
  A1[5,5]<-sum(T2plus*T2plus)
  A1[4,5]<-A1[5,4]<-sum(T2plus*T21plus)
  A1[4,4]<-sum(T12plus*T12plus)+sum(T21plus*T21plus)

  ## A3 is filled:
  A3[3,3]<-sum(T1)
  A3[5,5]<-sum(T2)
  A3[4,4]<-sum(T12)+sum(T21)
   

  ## A2 is filled:
  gamInv<-vecg^(-1)-1
  gi1<-gamInv[1];gi12<-gamInv[2];gi2<-gamInv[3]
  A2[1,1]<-sum(T1)*gi1
  A2[1,2]<-A2[2,1]<-sum(T12)*gi12
  A2[1,3]<-A2[3,1]<-sum(T1*(T1plus-1))*gi1
  A2[1,5]<-A2[5,1]<-sum(T21*T2plus)*gi12
  A2[1,4]<-A2[4,1]<-gi1*sum(T1*(T12plus))+gi12*sum(T21*(T21plus-1))
  A2[2,2]<-sum(T2)*gi2
  A2[2,3]<-A2[3,2]<-sum(T12*T1plus)*gi12
  A2[2,5]<-A2[5,2]<-sum(T2*(T2plus-1))*gi2
  A2[2,4]<-A2[4,2]<-gi2*sum(T2*(T21plus))+gi12*sum(T12*(T12plus-1))

  A2[3,3]<-gi1*sum(I1*outer(T1plus-1,T1plus-1))
  
  A2[3,5]<-A2[5,3]<- gi12*sum(I12*outer(T1plus,T2plus))
  A2[3,4]<-A2[4,3]<-gi1*sum(I1*outer(T1plus-1,T12plus))+gi12*sum(I12*outer(T1plus,T21plus-1))
  
  A2[5,5]<-gi2*sum(I2*outer(T2plus-1,T2plus-1))
  A2[4,5]<-A2[5,4]<-gi2*sum(I2*outer(T2plus-1,T21plus))+gi12*sum(I21*outer(T2plus,T12plus-1))
  
  A2[4,4]<-gi1*sum(I1*outer(T12plus,T12plus))+gi2*sum(I2*outer(T21plus,T21plus))+ gi12*sum(I12*outer(T12plus-1,T21plus-1))+gi12*sum(I21*outer(T21plus-1,T12plus-1))
  
  Sigma<-A1+A2+A3
  nam <- c(levels(marks(Xplus)), names(vecg))
  dnam <- list(nam, nam)
  
  if(spill) {
    # return internal data (with matrices unnormalised)
    dimnames(A1) <- dimnames(A2) <- dimnames(A3) <- dimnames(Sigma) <- dnam
    internals <- list(A1=A1, A2=A2, A3=A3, Sigma=Sigma)
    if(!spill.vc) return(internals)
  }
           
  ## Calculate variance-covariance
  # Normalise:
  A1    <- A1/areaW
  Sigma <- Sigma/areaW
  U <- checksolve(A1, matrix.action, , "variance")
  mat <- if(is.null(U)) matrix(NA, length(nam), length(nam)) else U%*%Sigma%*%U / areaW
  dimnames(mat) <- dnam

  if(spill.vc) return(list(varcov=mat, internals=internals))
  return(mat)
}

# Convert the first p rows & columns of variance matrix x
# to variances of treatment contrasts
contrastmatrix <- function(x,p){
  mat <- x
  ## Correct column and row 1:
  for(i in 2:p){
    mat[1,i] <- mat[i,1] <- x[1,i]-x[1,1]
  }
  ## Correct columns and rows 2,...,p:
  for(i in 2:p){
    for(j in 2:p){
      mat[i,j] <- x[1,1]-x[1,i]-x[1,j]+x[i,j]
    }
    for(j in (p+1):ncol(x)){
      mat[i,j] <- mat[j,i] <- x[i,j]-x[1,j]
    }
  }
  mat
}


checksolve <- function(M, action, descrip, target="") {
  Mname <- short.deparse(substitute(M))
  Minv <- try(solve(M), silent=(action=="silent"))
  if(!inherits(Minv, "try-error"))
    return(Minv)
  if(missing(descrip))
    descrip <- paste("the matrix", sQuote(Mname))
  whinge <- paste0("Cannot compute ", target, ": ", descrip, " is singular")
  switch(action,
         fatal=stop(whinge, call.=FALSE),
         warn= warning(whinge, call.=FALSE),
         silent={})
  return(NULL)
}

vcov.ppm
}
)

suffloc <- function(object) {
  verifyclass(object, "ppm")
  if(!is.poisson(object))
    stop("Internals not available for Gibbs models")
  return(vcov(object, what="internals")$suff)
}
