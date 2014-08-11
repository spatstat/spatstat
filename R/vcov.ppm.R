#
# Asymptotic covariance & correlation matrices
# and Fisher information matrix
# for ppm objects
#
#  $Revision: 1.78 $  $Date: 2013/05/23 08:05:06 $
#

vcov.ppm <- local({

vcov.ppm <- function(object, ..., what="vcov", verbose=TRUE,
                     gam.action=c("warn", "fatal", "silent"),
                     matrix.action=c("warn", "fatal", "silent"),
                     hessian=FALSE) {
  verifyclass(object, "ppm")
  argh <- list(...)
  
  gam.action <- match.arg(gam.action)
  matrix.action <- match.arg(matrix.action)

  stopifnot(length(what) == 1 && is.character(what))
  what.options <- c("vcov", "corr", "fisher", "Fisher", "internals")
  what.map     <- c("vcov", "corr", "fisher", "fisher", "internals")
  if(is.na(m <- pmatch(what, what.options)))
    stop(paste("Unrecognised option: what=", sQuote(what)))
  what <- what.map[m]

  # nonstandard calculations (hack) 
  generic.triggers <- c("A1", "A1dummy", "new.coef", "matwt", "saveterms",
                        "parallel")
  nonstandard <- any(generic.triggers %in% names(argh))
  saveterms <- identical(resolve.1.default("saveterms", argh), TRUE)
  
  # Fisher information *may* be contained in object
  fisher <- object$fisher
  varcov <- object$varcov
  
  # Do we need to go into the guts?
  needguts <-
    (is.null(fisher) && what=="fisher") ||
    (is.null(varcov) && what %in% c("vcov", "corr")) ||
    (what == "internals") || nonstandard

  # In general it is not true that varcov = solve(fisher)
  # because we might use different estimators,
  # or the parameters might be a subset of the canonical parameter

  if(needguts) {
    # warn if fitted model was obtained using GAM
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
    # ++++ perform main calculation ++++
    if(is.poisson(object) || hessian) {
      # Poisson model, or Hessian of Gibbs model
      results <- vcalcPois(object, ..., what=what,
                           matrix.action=matrix.action,
                           verbose=verbose, fisher=fisher)
    } else {
      # Gibbs model 
      results <- vcalcGibbs(object, ..., what=what,
                            matrix.action=matrix.action)
    }
    if(is.null(results))
      return(NULL)
    varcov <- results$varcov
    fisher <- results$fisher
    other  <- results$other
  }
  
  if(what %in% c("vcov", "corr") && is.null(varcov)) {
    # Need variance-covariance matrix.
    # Derive from Fisher information
    varcov <- checksolve(fisher, matrix.action,
                         "Fisher information matrix",
                         "variance")
    if(is.null(varcov)) return(NULL)
  }
         
  switch(what,
         fisher = { return(fisher) },
         vcov   = { return(varcov) },
         corr={
           sd <- sqrt(diag(varcov))
           return(varcov / outer(sd, sd, "*"))
         },
         internals = {
           return(append(list(fisher=fisher), other))
         })
}

# ................  variance calculation for Poisson models  .............

vcalcPois <- function(object, ...,
                      what = c("vcov", "corr", "fisher", "internals"),
                      matrix.action=c("warn", "fatal", "silent"),
                      method=c("C", "interpreted"),
                      verbose=TRUE,
                      fisher=NULL, 
                      matwt=NULL, new.coef=NULL,
                      saveterms=FALSE) {
  # variance-covariance matrix of Poisson model,
  # or Hessian of Gibbs model
  what <- match.arg(what)
  method <- match.arg(method)
  matrix.action <- match.arg(matrix.action)
  nonstandard <- !is.null(matwt) || !is.null(new.coef) || saveterms
  # compute Fisher information (Poincare information) if not known
  if(is.null(fisher) || nonstandard) {
    gf <- getglmfit(object)
    # we need a glm or gam
    if(is.null(gf)) {
      if(verbose) 
        warning("Refitting the model using GLM/GAM")
      object <- update(object, forcefit=TRUE)
      gf <- getglmfit(object)
      if(is.null(gf))
        stop("Internal error - refitting did not yield a glm object")
    }
    # compute fitted intensity and sufficient statistic
    ltype <- if(is.poisson(object)) "trend" else "lambda"
    lambda <- fitted(object, type=ltype, new.coef=new.coef, check=FALSE)
    mom <- model.matrix(object)
    nmom <- nrow(mom)
    Q <- quad.ppm(object)
    wt <- w.quad(Q)
    ok <- getglmsubset(object)
    Z  <- is.data(Q)
    # save them
    if(what == "internals") {
      internals <-
        if(!saveterms) list(suff=mom) else
      list(suff=mom, mom=mom, lambda=lambda, Z=Z, ok=ok)
    }
    # Now restrict all terms to the domain of the pseudolikelihood
    lambda <- lambda[ok]
    mom <- mom[ok, , drop=FALSE]
    wt <- wt[ok]
    Z <- Z[ok]
    # apply weights to rows of model matrix - temporary hack
    if(!is.null(matwt)) {
      nwt <- length(matwt)
      if(nwt == nmom) {
        # matwt matches original quadrature scheme - trim it
        matwt <- matwt[ok]
      } else if(nwt != sum(ok))
        stop("Hack argument matwt has incompatible length")
      mom <- matwt * mom
    }
    # compute Fisher information
    switch(method,
           C = {
             fisher <- sumouter(mom, lambda * wt)
           },
           interpreted = {
             for(i in 1:nrow(mom)) {
               ro <- mom[i, ]
               v <- outer(ro, ro, "*") * lambda[i] * wt[i]
               if(!any(is.na(v)))
                 fisher <- fisher + v
             }
             momnames <- dimnames(mom)[[2]]
             dimnames(fisher) <- list(momnames, momnames)
           })
  }
  switch(what,
         fisher = {
           result <- list(fisher=fisher)
         },
         corr = ,
         vcov = {
           # Derive variance-covariance from Fisher info
           varcov <- checksolve(fisher, matrix.action,
                                "Fisher information matrix",
                                "variance")
           result <- list(fisher=fisher, varcov=varcov)
         },
         internals = {
           result <- list(fisher=fisher, other=internals)
         })
  return(result)
}


# ...................... vcov calculation for Gibbs models ....................

vcalcGibbs <- function(fit, ...,
                       what = c("vcov", "corr", "fisher", "internals"),
                       generic=FALSE) {
  verifyclass(fit, "ppm")
  what <- match.arg(what)
  # decide whether to use the generic, slower algorithm
  generic.triggers <- c("A1", "A1dummy", "new.coef", "matwt", "saveterms",
                        "parallel")
  use.generic <-
    generic ||
  !is.stationary(fit) ||
  (fit$method == "logi" && ("marks" %in% variablesinformula(fit$trend))) ||
  (fit$method != "logi" && has.offset(fit)) ||
  (fit$method == "logi" && has.offset.term(fit)) ||
  !(fit$correction == "border" && fit$rbord == reach(fit)) ||
  any(generic.triggers %in% names(list(...)))
  !identical(options("contrasts")[[1]],
             c(unordered="contr.treatment",
               ordered="contr.poly"))
  #
  varcov <- if(use.generic) vcovGibbsAJB(fit, ...) else vcovGibbsEgeJF(fit, ...)
  if(is.null(varcov))
    return(NULL)
  # varcov is the variance-covariance matrix, with attributes
  Alist <- attr(varcov, "A")
  saved.terms <- attr(varcov, "saved.terms")
  internals <- append(Alist, saved.terms)
  # wipe attributes
  attr(varcov, "A") <- attr(varcov, "saved.terms") <- NULL
  #
   if(what %in% c("fisher", "internals")) {
    # compute fisher information from varcov
    fisher <- checksolve(varcov, matrix.action,
                         "variance-covariance matrix",
                         "Fisher information" )
   } else fisher <- NULL
  if(what == "internals") {
    # ensure sufficient statistic is calculated
    mom <- saved.terms$mom
    if(is.null(mom))
      mom <- model.matrix(fit)
    internals <- append(internals, list(suff=mom))
  } 
  return(list(varcov=varcov, fisher=fisher, other=internals))
}

## vcovGibbsAJB
## Adrian's code, modified by Ege to handle logistic case as well
##                modified by Adrian to handle 'new.coef'
## 2012/10/26, bugfixes by Ege for the logistic case

vcovGibbsAJB <- function(model,
                         ...,
                         matrix.action=c("warn", "fatal", "silent"),
                         algorithm=c("vectorclip", "vector", "basic"),
                         A1 = NULL,
                         A1dummy = FALSE,
                         matwt = NULL, new.coef = NULL,
                         saveterms = FALSE,
                         parallel = FALSE
                         ) {
  stopifnot(is.ppm(model))
  matrix.action <- match.arg(matrix.action)
  algorithm <- match.arg(algorithm)
  if(reweighting <- !is.null(matwt)) 
    stopifnot(is.numeric(matwt) && is.vector(matwt))
  logi <- model$method=="logi"
  old.coef <- coef(model)
  use.coef <- if(!is.null(new.coef)) new.coef else old.coef
  p <- length(old.coef)
  if(p == 0) {
    # this probably can't happen
    return(matrix(, 0, 0))
  }
  pnames <- names(old.coef)
  dnames <- list(pnames, pnames)
  saved.terms <- if(saveterms) list() else NULL
  #
  sumobj <- summary(model, quick="entries")
  correction <- model$correction
  rbord      <- model$rbord
  Q <- quad.ppm(model)
  D <- dummy.ppm(model)
  rho <- model$internal$logistic$rho
  ## If dummy intensity rho is unknown we estimate it
  if(is.null(rho))
     rho <- npoints(D)/(area.owin(D)*markspace.integral(D))
  X <- data.ppm(model)
  Z <- is.data(Q)
  W <- as.owin(model)
  #
  # determine which quadrature points contributed to the
  # sum/integral in the pseudolikelihood
  # (e.g. some points may be excluded by the border correction)
  okall <- getglmsubset(model)
  # data only:
  ok <- okall[Z]
  nX <- npoints(X)
  # conditional intensity lambda(X[i] | X) = lambda(X[i] | X[-i])
  # data and dummy:
  lamall <- fitted(model, check = FALSE, new.coef = new.coef)
  # data only:
  lam <- lamall[Z]
  # sufficient statistic h(X[i] | X) = h(X[i] | X[-i])
  # data and dummy:
  mall <- model.matrix(model)
  # save
  if(saveterms) 
    saved.terms <- append(saved.terms,
                          list(mom=mall, lambda=lamall, Z=Z, ok=okall,
                               matwt=matwt))
  if(reweighting) {
    # each column of the model matrix is multiplied by 'matwt'
    check.nvector(matwt, nrow(mall), things="quadrature points")
    mall <- mall * matwt
  }
  mokall <- mall[okall, , drop=FALSE]
  # data only:
  m <- mall[Z, , drop=FALSE]
  mok <- m[ok, , drop=FALSE]
  matwtX <- matwt[Z]

  # ^^^^^^^^^^^^^^^^ First order (sensitivity) matrices A1, S
  
  # logistic 
  if(logi){
    # Sensitivity matrix S for logistic case
    Slog <- sumouter(mokall, w = lamall[okall]*rho/(lamall[okall]+rho)^2)
    dimnames(Slog) <- dnames
    # A1 matrix for logistic case
    A1log <- sumouter(mokall, w = lamall[okall]*rho*rho/(lamall[okall]+rho)^3)
    dimnames(A1log) <- dnames
  }
  # Sensitivity matrix for MPLE case (= A1)
  if(is.null(A1)) {
    if(A1dummy){
#    A1 <- sumouter(mokall, w = lamall[okall]/(lamall[okall]+rho))
      A1 <- sumouter(mokall, w = (lamall * w.quad(Q))[okall])
    } else{
      A1 <- sumouter(mok)
    }
  } else {
    stopifnot(is.matrix(A1))
    if(!all(dim(A1) == p))
      stop(paste("Matrix A1 has wrong dimensions:",
                 prange(dim(A1)), "!=", prange(c(p, p))))
  }
  dimnames(A1) <- dnames

  # ^^^^^^^^^^ Second order interaction effects A2, A3
  
  # ^^^^^^^^^^^^^^^^^^^^ `parallel' evaluation

  need.loop <- TRUE
  if(parallel && !logi && require(tensor)) {
    # compute second order difference
    #  ddS[i,j,] = h(X[i] | X) - h(X[i] | X[-j])
    ddS <- deltasuffstat(model, restrict=FALSE)
    if(!is.null(ddS)) {
      need.loop <- FALSE
      # rearrange so that
      #  ddS[ ,i,j] = h(X[i] | X) - h(X[i] | X[-j])
      ddS <- aperm(ddS, c(3,2,1))
      # now compute sum_{i,j} for i != j
      # outer(ddS[,i,j], ddS[,j,i])
      ddSok <- ddS[ , ok, ok, drop=FALSE]
      A3 <- sumsymouter(ddSok)
      # momdif[ ,i,j] = h(X[i] | X[-j])
      momdif <- array(t(m), dim=c(p, nX, nX)) - ddS
      # lamdif[i,j] = lambda(X[i] | X[-j])
      lamdif <- matrix(lam, nX, nX) * exp(tensor(-use.coef, ddS, 1, 1))
      #   pairweight[i,j] = lamdif[i,j]/lambda[i] - 1 
      pairweight <- lamdif / lam - 1
      # now compute sum_{i,j} for i != j
      # pairweight[i,j] * outer(momdif[,i,j], momdif[,j,i])
      # for data points that contributed to the pseudolikelihood
      momdifok <- momdif[ , ok, ok, drop=FALSE]
      A2 <- sumsymouter(momdifok, w=pairweight[ok, ok])
    } else warning("parallel option not available - reverting to loop")
  }
  
  # ^^^^^^^^^^^^^^^^^^^^ loop evaluation

  if(need.loop) {
    
    A2 <- A3 <- matrix(0, p, p, dimnames=dnames)
    if(logi)
      A2log <- A3log <- matrix(0, p, p, dimnames=dnames)
    
    if(saveterms) {
      # *initialise* matrices 
      #  lamdif[i,j] = lambda(X[i] | X[-j]) = lambda(X[i] | X[-c(i,j)])
      lamdif <- matrix(lam, nX, nX)
      #  momdif[ ,i,j] = h(X[i] | X[-j]) = h(X[i] | X[-c(i,j)])
      momdif <- array(t(m), dim=c(p, nX, nX))
    }
  
    # identify close pairs
    R <- reach(model, epsilon=1e-2)
    if(is.finite(R)) {
#    if(R == 0 && !logi) return(A1)    ... redundant
#    if(R == 0 && logi) return(A1log)  ... redundant
      cl <- closepairs(X, R, what="indices")
      I <- cl$i
      J <- cl$j
      if(algorithm == "vectorclip") {
        cl2 <- closepairs(X, 2*R, what="indices")
        I2 <- cl2$i
        J2 <- cl2$j
      }
    } else {
      # either infinite reach, or something wrong
      IJ <- expand.grid(I=1:nX, J=1:nX)
      IJ <- subset(IJ, I != J)
      I2 <- I <- IJ$I
      J2 <- J <- IJ$J
    }
    # filter:  I and J must both belong to the nominated subset 
    okIJ <- ok[I] & ok[J]
    I <- I[okIJ]
    J <- J[okIJ]
    #
    if(length(I) > 0 && length(J) > 0) {
      # .............. loop over pairs ........................
      # The following ensures that 'empty' and 'X' have compatible marks 
      empty <- X[integer(0)]
      # make an empty 'equalpairs' matrix
      nonE <- matrix(, nrow=0, ncol=2)
      # Run through pairs
      switch(algorithm,
             basic={
               for(i in unique(I)) {
                 Xi <- X[i]
                 Ji <- unique(J[I==i])
                 if((nJi <- length(Ji)) > 0) {
                   for(k in 1:nJi) {
                     j <- Ji[k]
                     X.ij <- X[-c(i,j)]
                     # compute conditional intensity
                     #    lambda(X[j] | X[-i]) = lambda(X[j] | X[-c(i,j)]
                     plamj.i <- predict(model, type="cif",
                                        locations=X[j], X=X.ij,
                                        check = FALSE,
                                        new.coef = new.coef,
                                        sumobj = sumobj, E=nonE)
                     # corresponding values of sufficient statistic 
                     #    h(X[j] | X[-i]) = h(X[j] | X[-c(i,j)]
                     pmj.i <- partialModelMatrix(X.ij, X[j], model)[nX-1, ]
                     # conditional intensity and sufficient statistic
                     # in reverse order
                     #    lambda(X[i] | X[-j]) = lambda(X[i] | X[-c(i,j)]
                     plami.j <- predict(model, type="cif",
                                        locations=X[i], X=X.ij,
                                        check = FALSE,
                                        new.coef = new.coef,
                                        sumobj = sumobj, E=nonE)
                     pmi.j <- partialModelMatrix(X.ij, Xi, model)[nX-1, ]
                     # 
                     if(reweighting) {
                       pmj.i <- pmj.i * matwtX[j]
                       pmi.j <- pmi.j * matwtX[i]
                     }
                     if(saveterms) {
                       lamdif[i,j] <- plami.j
                       momdif[ , i, j] <- pmi.j
                       lamdif[j,i] <- plamj.i
                       momdif[ , j, i] <- pmj.i
                     }
                     # increment A2, A3
                     wt <- plami.j / lam[i] - 1
                     A2 <- A2 + wt * outer(pmi.j, pmj.i)
                     if(logi)
                     A2log <- A2log +
                       wt * rho/(plami.j+rho) *
                         rho/(plamj.i+rho) * outer(pmi.j, pmj.i)
                   # delta sufficient statistic
                   # delta_i h(X[j] | X[-c(i,j)])
                   # = h(X[j] | X[-j]) - h(X[j] | X[-c(i,j)])
                   # = h(X[j] | X) - h(X[j] | X[-i])
                   # delta_j h(X[i] | X[-c(i,j)])
                   # = h(X[i] | X[-i]) - h(X[i] | X[-c(i,j)])
                   # = h(X[i] | X) - h(X[i] | X[-j])
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
               # --------- faster algorithm using vector functions ------------
               for(i in unique(I)) {
                 Ji <- unique(J[I==i])
                 nJi <- length(Ji)
                 if(nJi > 0) {
                   Xi <- X[i]
                   # neighbours of X[i]
                   XJi <- X[Ji]
                   # all points other than X[i]
                   X.i <- X[-i]
                   # index of XJi in X.i
                   J.i <- Ji - (Ji > i)
                   # equalpairs matrix
                   E.i <- cbind(J.i, seq_len(nJi))
                   # compute conditional intensity
                   #   lambda(X[j] | X[-i]) = lambda(X[j] | X[-c(i,j)]
                   # for all j
                   plamj <- predict(model, type="cif",
                                    locations=XJi, X=X.i,
                                    check = FALSE,
                                    new.coef = new.coef,
                                    sumobj=sumobj, E=E.i)
                   # corresponding values of sufficient statistic 
                   #    h(X[j] | X[-i]) = h(X[j] | X[-c(i,j)]
                   # for all j
                 pmj <- partialModelMatrix(X.i, empty, model)[J.i, , drop=FALSE]
                 #
                 # conditional intensity & sufficient statistic in reverse order
                 #    lambda(X[i] | X[-j]) = lambda(X[i] | X[-c(i,j)]
                 # for all j
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
                   #
                   if(reweighting) {
                     pmj <- pmj * matwtX[Ji]
                     pmi <- pmi * matwtX[i]
                   }
                   if(saveterms) {
                     lamdif[Ji, i] <- plamj
                     momdif[ , Ji, i] <- t(pmj)
                     lamdif[i,Ji] <- plami
                     momdif[ , i, Ji] <- t(pmi)
                   }
                   # increment A2, A3
                   wt <- plami / lam[i] - 1
                   for(k in 1:nJi) {
                     j <- Ji[k]
                     A2 <- A2 + wt[k] * outer(pmi[k,], pmj[k,])
                     if(logi)
                     A2log <- A2log + wt[k] * rho/(plami[k]+rho) *
                       rho/(plamj[k]+rho) * outer(pmi[k,], pmj[k,])
                   # delta sufficient statistic
                   # delta_i h(X[j] | X[-c(i,j)])
                   # = h(X[j] | X[-j]) - h(X[j] | X[-c(i,j)])
                   # = h(X[j] | X) - h(X[j] | X[-i])
                   # delta_j h(X[i] | X[-c(i,j)])
                   # = h(X[i] | X[-i]) - h(X[i] | X[-c(i,j)])
                   # = h(X[i] | X) - h(X[i] | X[-j])
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
               # --------- faster version of 'vector' algorithm
               # --------  by removing non-interacting points of X
               for(i in unique(I)) {
                 # all points within 2R
                 J2i <- unique(J2[I2==i])
                 # all points within R
                 Ji  <- unique(J[I==i])
                 nJi <- length(Ji)
                 if(nJi > 0) {
                   Xi <- X[i]
                   # neighbours of X[i]
                   XJi <- X[Ji]
                   # replace X[-i] by X[-i] \cap b(0, 2R)
                   X.i <- X[J2i]
                   nX.i <- length(J2i)
                   # index of XJi in X.i
                   J.i <- match(Ji, J2i)
                   if(any(is.na(J.i)))
                     stop("Internal error: Ji not a subset of J2i")
                   # equalpairs matrix
                   E.i <- cbind(J.i, seq_len(nJi))
                   # compute conditional intensity
                   #   lambda(X[j] | X[-i]) = lambda(X[j] | X[-c(i,j)]
                   # for all j
                   plamj <- predict(model, type="cif",
                                    locations=XJi, X=X.i,
                                    check = FALSE,
                                    new.coef = new.coef,
                                    sumobj = sumobj, E=E.i)
                   # corresponding values of sufficient statistic 
                   #    h(X[j] | X[-i]) = h(X[j] | X[-c(i,j)]
                   # for all j
                 pmj <- partialModelMatrix(X.i, empty, model)[J.i, , drop=FALSE]
                 #
                 # conditional intensity & sufficient statistic in reverse order
                 #    lambda(X[i] | X[-j]) = lambda(X[i] | X[-c(i,j)]
                 # for all j
                   plami <- numeric(nJi)
                   pmi <- matrix(, nJi, p)
                   for(k in 1:nJi) {
                     j <- Ji[k]
                     # X.ij <- X[-c(i,j)]
                     X.ij <- X.i[-J.i[k]]
                     plami[k] <- predict(model, type="cif",
                                         locations=Xi, X=X.ij,
                                         check = FALSE,
                                         new.coef = new.coef,
                                         sumobj = sumobj, E=nonE)
                     pmi[k, ] <- partialModelMatrix(X.ij, Xi, model)[nX.i, ]
                   }
                   #
                   if(reweighting) {
                     pmj <- pmj * matwtX[Ji]
                     pmi <- pmi * matwtX[i]
                   }
                   if(saveterms) {
                     lamdif[Ji, i] <- plamj
                     momdif[ , Ji, i] <- t(pmj)
                     lamdif[i,Ji] <- plami
                     momdif[ , i, Ji] <- t(pmi)
                   }
                   # increment A2, A3
                   wt <- plami / lam[i] - 1
                   for(k in 1:nJi) {
                     j <- Ji[k]
                     A2 <- A2 + wt[k] * outer(pmi[k,], pmj[k,])
                     if(logi)
                       A2log <- A2log + wt[k] * rho/(plami[k]+rho) *
                         rho/(plamj[k]+rho) * outer(pmi[k,], pmj[k,])
                   # delta sufficient statistic
                   # delta_i h(X[j] | X[-c(i,j)])
                   # = h(X[j] | X[-j]) - h(X[j] | X[-c(i,j)])
                   # = h(X[j] | X) - h(X[j] | X[-i])
                   # delta_j h(X[i] | X[-c(i,j)])
                   # = h(X[i] | X[-i]) - h(X[i] | X[-c(i,j)])
                   # = h(X[i] | X) - h(X[i] | X[-j])
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
  # ......... end of loop computation ...............

  if(saveterms) 
    saved.terms <- append(saved.terms, list(lamdif=lamdif, momdif=momdif))
  # Normalise by eroded area
  areaW <- if(correction == "border") eroded.areas(W, rbord) else area.owin(W)
  #
  A1 <- A1/areaW
  A2 <- A2/areaW
  A3 <- A3/areaW
  if(logi){
    A1log <- A1log/areaW
    A2log <- A2log/areaW
    A3log <- A3log/areaW
    Slog <- Slog/areaW
  }
  ## Matrix Sigma (A1+A2+A3):
  Sigma <- A1+A2+A3
  # Enforce exact symmetry of A1 and Sigma
  A1 <- (A1 + t(A1))/2
  Sigma <- (Sigma + t(Sigma))/2
  ## Finally the result is calculated:
  U <- checksolve(A1, matrix.action, , "variance")
  if(is.null(U)) {
    mat <- matrix(NA, p, p)
  } else {
    mat <- U %*% Sigma %*% U / areaW
    dimnames(mat) <- dnames
  }
  # Add extra information
  Alist <- list(A1=A1, A2=A2, A3=A3, areaW=areaW, Sigma=Sigma)
  if(!logi) {
    attr(mat, "A") <- Alist    
    attr(mat, "saved.terms") <- saved.terms
    return(mat)
  }
  ###### Everything below is only computed for logistic fits #######
  
  Alist <- append(Alist,
                  list(A1log=A1log, A2log=A2log, A3log=A3log, Slog=Slog))

  ## Matrix Sigma1log (A1log+A2log+A3log):
  Sigma1log <- A1log+A2log+A3log
  ## Matrix Sigma2log (depends on dummy process type)
  switch(model$internal$logistic$how,
         poisson={
           Sigma2log <- sumouter(mokall, w = lamall[okall]*lamall[okall]*rho/(lamall[okall]+rho)^3)
           Sigma2log <- Sigma2log/areaW
         },
         binomial={
           Sigma2log <- sumouter(mokall, w = lamall[okall]*lamall[okall]*rho/(lamall[okall]+rho)^3)/areaW
           A1vec <- t(mokall) %*% (rho*lamall[okall]/(lamall[okall]+rho)^2/areaW)
           Sigma2log <- Sigma2log - A1vec%*%t(A1vec)/rho*areaW/sum(1/(lamall[okall]+rho))
         },
         stratrand={
           ### Dirty way of refitting model with new dummy pattern (should probably be done using call, eval, envir, etc.):
           D2 <- logi.dummy(X = X, type = "stratrand", nd = model$internal$logistic$args)
           Q2 <- quad(data=X, dummy=D2)
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
# Excised by AJB           
#           if(reweighting) {
#             mall2 <- mall2 * matwt
#           }
           okall2 <- getglmsubset(model2)

           # index vectors of stratrand cell indices of dummy points 
           inD <- model$internal$logistic$inD
           inD2 <- model2$internal$logistic$inD

           # Dummy points inside eroded window (for border correction)
           if(is.finite(R) && (correction == "border")){
             ii <- (bdist.points(D) >= R)
             ii2 <- (bdist.points(D2) >= R)
#  AJB:  faster and more consistent than:             
#             ii <- inside.owin(D, w = Wminus)
#             ii2 <- inside.owin(D2, w = Wminus)
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
           Sigma2log <- t(wlam-wlam2)%*%(wlam-wlam2)/(2*rho*rho*areaW)
         },
         stop("sorry - unrecognized dummy process in logistic fit")
         )
  ## Attaching to Sigma2log calculated above
  dimnames(Sigma2log) <- dnames

  ## Finally the result is calculated:
  Ulog <- checksolve(Slog, matrix.action, , "variance")
  if(is.null(Ulog)) {
    matlog <- matrix(NA, p, p)
  } else {
    matlog <- Ulog %*% (Sigma1log+Sigma2log) %*% Ulog / areaW
    dimnames(matlog) <- dnames
  }
  #
  attr(matlog, "A") <-
    append(Alist,
           list(Sigma1log=Sigma1log, Sigma2log=Sigma2log, mple=mat))
  attr(matlog, "saved.terms") <- saved.terms

  return(matlog)
}

# vcovGibbs from Ege Rubak and J-F Coeurjolly
## 2012/10/26, modified by Ege to handle logistic case as well

vcovGibbsEgeJF <- function(fit, ...,
                           special.alg = TRUE,
                           matrix.action=c("warn", "fatal", "silent")) {
  verifyclass(fit, "ppm")
  matrix.action <- match.arg(matrix.action)
  
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
                               matrix.action))
      },
           
      "Piecewise constant pairwise interaction process"={
        ## Only implemented for non-marked case:
        if(!marx)
          return(vcovPairPiece(Xplus,
                               fit$interaction$par$r,
                               exp(coef(fit)[-1]),
                               matrix.action))
      },

      "Multitype Strauss process"={
	matR <- fit$interaction$par$radii
        R <- c(matR[1,1], matR[1,2], matR[2,2])
        ## Only implemented for 2 types with equal interaction range:
        if(ncol(matR)==2 && marx){
          n <- length(theta)
          res <- vcovMultiStrauss(Xplus, R, exp(theta[c(n-2,n-1,n)]),
                                  matrix.action)
          res <- contrastmatrix(res, 2)
          dimnames(res) <- list(names(theta), names(theta))
          return(res)
        }
      }
    )
  }
  
  ## Matrix specifying equal points in the two patterns in the call to eval below:
  E <- matrix(rep.int(1:n, 2), ncol = 2)

  ## Eval. the interaction potential difference at all points (internal spatstat function):
  V1 <- fit$interaction$family$eval(Xplus, Xplus, E, fit$interaction$pot, fit$interaction$par, fit$correction)

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
    if(fit$interaction$family$name=="pairwise"){
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
  areaW <- area.owin(W)

  ## Interior points determined by bdist.points:
  IntPoints <- bdist.points(Xplus)>=R  
  X <- Xplus[IntPoints]
  
  ## Making a logical matrix, I, indicating R-close pairs which are in the interior:
  D <- pairdist(Xplus)
  diag(D) <- Inf
  I <- D<=R * outer(IntPoints,IntPoints)
  
  ## Matrix A1:
  A1 <- t(V1[IntPoints,])%*%V1[IntPoints,]/areaW

  ## Matrix A2:
  A2 <- matrix(0,p,p)
  for(k in 1:p){
    for(l in k:p){
      A2[k,l] <- A2[l,k] <- sum(I*V2[,,k]*frac*t(V2[,,l]))
    }
  }
  A2 <- A2/areaW
  
  ## Matrix A3:
  A3 <- matrix(0,p,p)
  for(k in 1:p){
    for(l in k:p){
      A3[k,l] <- A3[l,k] <- sum(I*dV[,,k]*t(dV[,,l]))
    }
  }
  A3 <- A3/areaW
  
  ## Matrix Sigma (A1+A2+A3):
  Sigma<-A1+A2+A3
  
  ## Finally the result is calculated:
  U <- checksolve(A1, matrix.action, , "variance")
  if(is.null(U)) return(NULL)

  mat <- U%*%Sigma%*%U / areaW

  ## Convert to treatment contrasts
 if(marx){
    mat <- contrastmatrix(mat, p1)
    dimnames(mat) <- list(names(theta), names(theta))
  }
  # save A1, A2, A3 
  dimnames(A1) <- dimnames(A2) <-
    dimnames(A3) <- list(names(theta), names(theta))
  attr(mat, "A") <- list(A1=A1, A2=A2, A3=A3)
  #
  ## Return result for standard ppm method:
  if(fit$method!="logi")
    return(mat)

  ########################################################################
  ###### The remainder is only executed when the method is logistic ######
  ########################################################################

  ### Most of this is copy/pasted from vcovGibbsAJB
  correction <- fit$correction
  Q <- quad.ppm(fit)
  D <- dummy.ppm(fit)
  rho <- fit$internal$logistic$rho
  ## If dummy intensity rho is unknown we estimate it
  if(is.null(rho))
     rho <- npoints(D)/(area.owin(D)*markspace.integral(D))
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
  Slog <- sumouter(mokall, w = lamall[okall]*rho/(lamall[okall]+rho)^2)/areaW
  A1log <- sumouter(mokall, w = lamall[okall]*rho*rho/(lamall[okall]+rho)^3)/areaW

  ## Define W1, W2 and dW for the logistic method based on V1, V2 and dV (frac is unchanged)
  lambda1 <- exp(rowSums(matrix(theta,n,p,byrow=TRUE)*V1))
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
  A2log <- A2log/areaW
  A3log <- A3log/areaW
  
  ## First variance component Sigma1log (A1log+A2log+A3log):
  Sigma1log <- A1log+A2log+A3log

  ## Matrix Sigma2log (depends on dummy process type)
  switch(fit$internal$logistic$how,
         poisson={
           Sigma2log <- sumouter(mokall, w = lamall[okall]*lamall[okall]*rho/(lamall[okall]+rho)^3)/areaW
         },
         binomial={
           Sigma2log <- sumouter(mokall, w = lamall[okall]*lamall[okall]*rho/(lamall[okall]+rho)^3)/areaW
           A1vec <- t(mokall) %*% (rho*lamall[okall]/(lamall[okall]+rho)^2/areaW)
           Sigma2log <- Sigma2log - A1vec%*%t(A1vec)/rho*areaW/sum(1/(lamall[okall]+rho))
         },
         stratrand={
           ### Dirty way of refitting model with new dummy pattern (should probably be done using call, eval, envir, etc.):
           D2 <- logi.dummy(X = X, type = "stratrand", nd = fit$internal$logistic$args)
           Q2 <- quad(data=X, dummy=D2)
           Q2$dummy$Dinfo <- D2$Dinfo
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
           Sigma2log <- t(wlam-wlam2)%*%(wlam-wlam2)/(2*rho*rho*areaW)
         },
         stop("sorry - unrecognized dummy process in logistic fit")
         )

  ## Finally the result is calculated:
  Ulog <- checksolve(Slog, matrix.action, , "variance")
  if(is.null(Ulog)) return(NULL)
  matlog <- Ulog %*% (Sigma1log+Sigma2log) %*% Ulog / areaW
  
  ## Attaching dimnames to all matrices
  dimnames(Sigma2log) <- dimnames(matlog) <- dimnames(Slog) <- dimnames(Sigma1log) <- dimnames(A1log) <- dimnames(A2log) <- dimnames(A3log) <- list(names(theta),names(theta))

  #
  attr(matlog, "A") <- list(A1log=A1log, A2log=A2log, A3log=A3log, Slog=Slog, Sigma1log=Sigma1log, Sigma2log=Sigma2log, mple=mat)

  return(matlog)
}

vcovPairPiece <- function(Xplus, R, gam, matrix.action){
  ## R is  the  vector of breaks (R[length(R)]= range of the pp.
  ## gam is the vector of weights
  Rmax <- R[length(R)]
  
  ## Xplus : point process observed in W+R
  ## Extracting the window and calculating area:
  Wplus<-as.owin(Xplus)
  W<-erosion.owin(Wplus,Rmax)
  areaW <- area.owin(W)

  ## Interior points determined by bdist.points:
  IntPoints <- bdist.points(Xplus)>=Rmax
  X <- Xplus[IntPoints]
  
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
     Tplus[,i]<-colSums(Iplus[[i]])[IntPoints]
     T[,i] <- colSums(I[[i]])
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
    A2[1,1]<-A2[1,1]+(gam[j-1]^(-1)-1)*sum(T[,j-1])
    for (l in (2:(p+1))){
      if (l==j) vj<-Tplus[,j-1]-1 else vj<-Tplus[,j-1]
	A2[1,j]<-A2[1,j]+(gam[l-1]^(-1)-1)*sum(T[,l-1]*(vj) )
    }
    A2[j,1]<-A2[1,j]
    for (k in (2:(p+1))){
      for (l in (2:(p+1))){
	if (l==j) vj<-Tplus[,j-1]-1 else vj<-Tplus[,j-1]
	if (l==k) vk<-Tplus[,k-1]-1 else vk<-Tplus[,k-1]

	A2[j,k]<-A2[j,k]+ (gam[l-1]^(-1)-1)*sum(I[[l-1]]*outer(vj,vk))
      }
    }

  }
  A1<-A1/areaW
  A2<-A2/areaW
  A3<-A3/areaW
  ## browser()
  Sigma<-A1+A2+A3
  ## Finally the result is calculated:
  U <- checksolve(A1, matrix.action, , "variance")
  if(is.null(U)) return(NULL)

  mat <- U%*%Sigma%*%U / areaW
  nam <- c("(Intercept)", names(gam))
  dimnames(mat) <- list(nam, nam)
  dimnames(A1) <- dimnames(A2) <- dimnames(A3) <- list(nam, nam)
  attr(mat, "A") <- list(A1=A1, A2=A2, A3=A3)
  return(mat)
}

vcovMultiStrauss <- function(Xplus, vecR, vecg, matrix.action){
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
  areaW <- area.owin(W)
  X1plus<-Xplus[Xplus$marks==levels(Xplus$marks)[1]]
  X2plus<-Xplus[Xplus$marks==levels(Xplus$marks)[2]]

  ## Interior points determined by bdist.points:
  IntPoints1 <- bdist.points(X1plus)>=R
  IntPoints2 <- bdist.points(X2plus)>=R
  X1 <- X1plus[IntPoints1]
  X2 <- X2plus[IntPoints2]

  ## Matrix D with pairwise distances between points and infinite distance
  ## between a point and itself:

  D1plus<-pairdist(X1plus)
  D1 <- pairdist(X1)
  diag(D1) <- diag(D1plus) <- Inf
  
  D2plus<-pairdist(X2plus)
  D2 <- pairdist(X2)
  diag(D2) <- diag(D2plus) <- Inf
  
  D12plus<-crossdist(X1,X2plus)  
  T12plus<-rowSums(D12plus<=R12)
  D21plus<-crossdist(X2,X1plus) 
  T21plus<-rowSums(D21plus<=R12)
  
  I12<-crossdist(X1,X2)<=R12
  I21<-crossdist(X2,X1)<=R12
  T12<-rowSums( I12)  
  T21<-rowSums(I21)
  ## logical matrix, I, indicating R-close pairs:
  I1plus<- D1plus <=R11
  I1 <- D1<=R11
  I2plus<- D2plus <=R22
  I2 <- D2<=R22
  ## Vector T with the number of $R$-close neighbours to each point:
  T1plus<-colSums(I1plus)[IntPoints1]
  T1 <- colSums(I1)
  T2plus<-colSums(I2plus)[IntPoints2]
  T2 <- colSums(I2)

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
  
  #browser()
  A1<-A1/areaW
  A2<-A2/areaW
  A3<-A3/areaW
  #browser()
  
  Sigma<-A1+A2+A3
  ## Finally the result is calculated:
  U <- checksolve(A1, matrix.action, , "variance")
  if(is.null(U)) return(NULL)

  mat <- U%*%Sigma%*%U / areaW

  nam <- c(levels(marks(Xplus)), names(vecg))
  dimnames(mat) <- list(nam, nam)
  dimnames(A1) <- dimnames(A2) <- dimnames(A3) <- list(nam, nam)
  
  attr(mat, "A") <- list(A1=A1, A2=A2, A3=A3)
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
    descrip <- paste("the matrix", Mname)
  whinge <- paste("Cannot compute",
                  paste(target, ":", sep=""),
                  descrip, "is singular")
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
