#
# kppm.R
#
# kluster/kox point process models
#
# $Revision: 1.88 $ $Date: 2014/08/27 09:39:08 $
#

kppm <- function(X, ...) {
  UseMethod("kppm")
}


kppm.formula <-
  function(X, clusters = c("Thomas","MatClust","Cauchy","VarGamma","LGCP"),
           ..., data=NULL) {
  ## remember call
  callstring <- short.deparse(sys.call())
  ## cl <- match.call()

  ########### INTERPRET FORMULA ##############################
  
  if(!inherits(X, "formula"))
    stop(paste("Argument 'X' should be a formula"))
  formula <- X
  
  if(spatstat.options("expand.polynom"))
    formula <- expand.polynom(formula)

  ## check formula has LHS and RHS. Extract them
  if(length(formula) < 3)
    stop(paste("Formula must have a left hand side"))
  Yexpr <- lhs <- formula[[2]]
  trend <- rhs <- formula[c(1,3)]
  
  ## FIT #######################################
  thecall <- call("kppm", X=Yexpr, trend=trend,
                  data=data, clusters=clusters)
  ncall <- length(thecall)
  argh <- list(...)
  nargh <- length(argh)
  if(nargh > 0) {
    thecall[ncall + 1:nargh] <- argh
    names(thecall)[ncall + 1:nargh] <- names(argh)
  }
  result <- eval(thecall, parent.frame())

  if(!("callstring" %in% names(list(...))))
    result$callstring <- callstring
  
  return(result)
}

kppm.ppp <- kppm.quad <-
  function(X, trend = ~1,
           clusters = c("Thomas","MatClust","Cauchy","VarGamma","LGCP"),
           data=NULL,
           ...,
           covariates = data,
           method = c("mincon", "clik2", "palm"),
           improve.type = c("none", "clik1", "wclik1", "quasi"),
           improve.args = list(),
           weightfun=NULL,
           control=list(),
           statistic="K",
           statargs=list(),
           rmax = NULL,
           covfunargs=NULL,
           use.gam=FALSE,
           nd=NULL, eps=NULL) {
  Xname <- short.deparse(substitute(X))
  clusters <- match.arg(clusters)
  improve.type <- match.arg(improve.type)
  method <- match.arg(method)
  if(method == "mincon")
    statistic <- pickoption("summary statistic", statistic,
                            c(K="K", g="pcf", pcf="pcf"))
  isquad <- inherits(X, "quad")
  if(!is.ppp(X) && !isquad)
    stop("X should be a point pattern (ppp) or quadrature scheme (quad)")
  if(is.marked(X))
    stop("Sorry, cannot handle marked point patterns")
  po <- ppm(Q=X, trend=trend, covariates=covariates,
            forcefit=TRUE, rename.intercept=FALSE,
            covfunargs=covfunargs, use.gam=use.gam, nd=nd, eps=eps)
  XX <- if(isquad) X$data else X
  # fit
  out <- switch(method,
         mincon = kppmMinCon(X=XX, Xname=Xname, po=po, clusters=clusters,
                             statistic=statistic, statargs=statargs,
                             control=control, rmax=rmax, ...),
         clik2   = kppmComLik(X=XX, Xname=Xname, po=po, clusters=clusters,
                             control=control, weightfun=weightfun, 
                             rmax=rmax, ...),
         palm   = kppmPalmLik(X=XX, Xname=Xname, po=po, clusters=clusters,
                             control=control, weightfun=weightfun, 
                             rmax=rmax, ...))
  #
  class(out) <- c("kppm", class(out))

  # Update intensity estimate with improve.kppm if necessary:
  if(improve.type != "none")
    out <- do.call(improve.kppm,
                   append(list(object = out, type = improve.type),
                          improve.args))
  return(out)
}

kppmMinCon <- function(X, Xname, po, clusters, statistic, statargs, ...) {
  # Minimum contrast fit
  stationary <- is.stationary(po)
  # compute summary function
  if(stationary) {
    StatFun <- if(statistic == "K") "Kest" else "pcf"
    StatName <-
      if(statistic == "K") "K-function" else "pair correlation function"
    Stat <- do.call(StatFun,
                    resolve.defaults(list(X=X),
                                     statargs,
                                     list(correction="best")))
    lambda <- summary(po)$trend$value
  } else {
    StatFun <- if(statistic == "K") "Kinhom" else "pcfinhom"
    StatName <- if(statistic == "K") "inhomogeneous K-function" else
    "inhomogeneous pair correlation function"
    # compute intensity at high resolution if available
    w <- as.owin(po, from="covariates")
    if(!is.mask(w)) w <- NULL
    lambda <- predict(po, locations=w)
    Stat <- do.call(StatFun,
                    resolve.defaults(list(X=X, lambda=lambda),
                                     statargs,
                                     list(correction="best")))
  }
  # determine initial values of parameters
  selfstart <- spatstatClusterModelInfo(clusters)$selfstart
  startpar0 <- selfstart(X) 
  # fit
  switch(clusters,
         Thomas={
           FitFun <- if(statistic == "K") "thomas.estK" else "thomas.estpcf"
           mcfit <-
             do.call(FitFun,
                     resolve.defaults(
                                      list(X=Stat, lambda=lambda),
                                      list(...),
                                      list(startpar=startpar0, dataname=Xname)))
           # kappa = intensity of parents
           kappa <- mcfit$par[["kappa"]]
           # mu = mean number of points per cluster
           mu <- if(stationary) lambda/kappa else eval.im(lambda/kappa)
           isPCP <- TRUE
         },
         MatClust={
           FitFun <- if(statistic == "K") "matclust.estK" else "matclust.estpcf"
           mcfit <-
             do.call(FitFun,
                     resolve.defaults(
                                      list(X=Stat, lambda=lambda),
                                      list(...),
                                      list(startpar=startpar0, dataname=Xname)))
           # kappa = intensity of parents
           kappa <- mcfit$par[["kappa"]]
           # mu = mean number of points per cluster
           mu <- if(stationary) lambda/kappa else eval.im(lambda/kappa)
           isPCP <- TRUE
         },
         Cauchy = {
           FitFun <- if (statistic == "K") "cauchy.estK" else "cauchy.estpcf"
           mcfit <- do.call(FitFun,
                            resolve.defaults(list(X = Stat, 
                                                  lambda = lambda),
                                             list(...),
                                             list(startpar = startpar0,
                                                  dataname = Xname)))
           # kappa = intensity of parents
           kappa <- mcfit$par[["kappa"]]
           # mu = mean number of points per cluster
           mu <- if (stationary) lambda/kappa else eval.im(lambda/kappa)
           isPCP <- TRUE
         }, VarGamma = {
           FitFun <- if (statistic == "K") "vargamma.estK" else "vargamma.estpcf"
           mcfit <- do.call(FitFun,
                            resolve.defaults(list(X = Stat, 
                                                  lambda = lambda),
                                             list(...),
                                             list(startpar = startpar0, 
                                                  dataname = Xname, nu = 0.5)))
           kappa <- mcfit$par[["kappa"]]
           mu <- if (stationary) lambda/kappa else eval.im(lambda/kappa)
           isPCP <- TRUE
         },
         LGCP={
           FitFun <- if(statistic == "K") "lgcp.estK" else "lgcp.estpcf"
           mcfit <-
             do.call(FitFun,
                     resolve.defaults(
                                      list(X=Stat, lambda=lambda),
                                      list(...),
                                      list(startpar=startpar0, dataname=Xname)))
           sigma2 <- mcfit$par[["sigma2"]]
           # mu = mean of log intensity 
           mu <- if(stationary) log(lambda) - sigma2/2 else eval.im(log(lambda) - sigma2/2)
           isPCP <- FALSE
         })

  # all info that depends on the fitting method:
  Fit <- list(method    = "mincon",
              statistic = statistic,
              Stat      = Stat,
              StatFun   = StatFun,
              StatName  = StatName,
              FitFun    = FitFun,
              statargs  = statargs,
              mcfit     = mcfit)
  # results
  out <- list(Xname      = Xname,
              X          = X,
              stationary = stationary,
              clusters   = clusters,
              modelname  = mcfit$info$modelname,
              isPCP      = isPCP,
              po         = po,
              lambda     = lambda,
              mu         = mu,
              par        = mcfit$par,
              modelpar   = mcfit$modelpar,
              covmodel   = mcfit$covmodel,
              Fit        = Fit)
  return(out)
}

kppmComLik <- function(X, Xname, po, clusters, control, weightfun, rmax, ...) {
  W <- as.owin(X)
  if(is.null(rmax))
    rmax <- rmax.rule("K", W, intensity(X))
  # identify pairs of points that contribute
  cl <- closepairs(X, rmax)
  I <- cl$i
  J <- cl$j
  dIJ <- cl$d
  # compute weights for pairs of points
  if(is.function(weightfun)) {
    wIJ <- weightfun(dIJ)
    sumweight <- sum(wIJ)
  } else {
    npairs <- length(dIJ)
    wIJ <- rep.int(1, npairs)
    sumweight <- npairs
  }
  # convert window to mask, saving other arguments for later
  dcm <- do.call.matched("as.mask",
                         append(list(w=W), list(...)),
                         sieve=TRUE)
  M         <- dcm$result
  otherargs <- dcm$otherargs
  # compute intensity at pairs of data points
  # and c.d.f. of interpoint distance in window
  if(stationary <- is.stationary(po)) {
    # stationary unmarked Poisson process
    lambda <- intensity(X)
    lambdaIJ <- lambda^2
    # compute cdf of distance between two uniform random points in W
    g <- distcdf(W)
    # scaling constant is (area * intensity)^2
    gscale <- npoints(X)^2  
  } else {
    # compute fitted intensity at data points and in window
    lambdaX <- fitted(po, dataonly=TRUE)
    lambda <- lambdaM <- predict(po, locations=M)
    # lambda(x_i) * lambda(x_j)
    lambdaIJ <- lambdaX[I] * lambdaX[J]
    # compute cdf of distance between two random points in W
    # with density proportional to intensity function
    g <- distcdf(M, dW=lambdaM)
    # scaling constant is (integral of intensity)^2
    gscale <- integral.im(lambdaM)^2
  }
  # trim 'g' to [0, rmax] 
  g <- g[with(g, .x) <= rmax,]
  # get pair correlation function (etc) for model
  info <- spatstatClusterModelInfo(clusters)
  pcfun      <- info$pcf
  funaux     <- info$funaux
  selfstart  <- info$selfstart
  isPCP      <- info$isPCP
  parhandler <- info$parhandler
  modelname  <- info$modelname
  # Assemble information required for computing pair correlation
  pcfunargs <- list(funaux=funaux)
  if(is.function(parhandler)) {
    # Additional parameters of cluster model are required.
    # These may be given as individual arguments,
    # or in a list called 'covmodel'
    clustargs <- if("covmodel" %in% names(otherargs))
                 otherargs[["covmodel"]] else otherargs
    clargs <- do.call(parhandler, clustargs)
    pcfunargs <- append(clargs, pcfunargs)
  } else clargs <- NULL
  # determine starting parameter values
  startpar <- selfstart(X)
  # create local function to evaluate pair correlation
  #  (with additional parameters 'pcfunargs' in its environment)
  paco <- function(d, par) {
    do.call(pcfun, append(list(par=par, rvals=d), pcfunargs))
  }
  # define objective function 
  if(!is.function(weightfun)) {
    # pack up necessary information
    objargs <- list(dIJ=dIJ, sumweight=sumweight, g=g,
                    envir=environment(paco))
    # define objective function (with 'paco' in its environment)
    # Note that this is 1/2 of the log composite likelihood,
    # minus the constant term 
    #       sum(log(lambdaIJ)) - npairs * log(gscale)
    obj <- function(par, objargs) {
      with(objargs,
           sum(log(paco(dIJ, par)))
           - sumweight * log(unlist(stieltjes(paco, g, par=par))),
           enclos=objargs$envir)
    }
  } else {
    # create local function to evaluate  pair correlation(d) * weight(d)
    #  (with additional parameters 'pcfunargs', 'weightfun' in its environment)
    force(weightfun)
    wpaco <- function(d, par) {
      y <- do.call(pcfun, append(list(par=par, rvals=d), pcfunargs))
      w <- weightfun(d)
      return(y * w)
    }
    # pack up necessary information
    objargs <- list(dIJ=dIJ, wIJ=wIJ, sumweight=sumweight, g=g,
                    envir=environment(wpaco))
    # define objective function (with 'paco', 'wpaco' in its environment)
    # Note that this is 1/2 of the log composite likelihood,
    # minus the constant term 
    #       sum(wIJ * log(lambdaIJ)) - sumweight * log(gscale)
    obj <- function(par, objargs) {
      with(objargs,
           sum(wIJ * log(paco(dIJ, par)))
           - sumweight * log(unlist(stieltjes(wpaco, g, par=par))),
           enclos=objargs$envir)
    }
  }    
  # optimize it
  ctrl <- resolve.defaults(list(fnscale=-1), control, list(trace=0))
  opt <- optim(startpar, obj, objargs=objargs, control=ctrl)
  # raise warning/error if something went wrong
  signalStatus(optimStatus(opt), errors.only=TRUE)
  # meaningful model parameters
  optpar <- opt$par
  modelpar <- info$interpret(optpar, lambda)
  # infer parameter 'mu'
  if(isPCP) {
    # Poisson cluster process: extract parent intensity kappa
    kappa <- optpar[["kappa"]]
    # mu = mean cluster size
    mu <- if(stationary) lambda/kappa else eval.im(lambda/kappa)
  } else {
    # LGCP: extract variance parameter sigma2
    sigma2 <- optpar[["sigma2"]]
    # mu = mean of log intensity 
    mu <- if(stationary) log(lambda) - sigma2/2 else
          eval.im(log(lambda) - sigma2/2)    
  }
  # all info that depends on the fitting method:
  Fit <- list(method    = "clik2",
              clfit     = opt,
              weightfun = weightfun,
              rmax      = rmax,
              objfun    = obj,
              objargs   = objargs)
  # pack up
  result <- list(Xname      = Xname,
                 X          = X,
                 stationary = stationary,
                 clusters   = clusters,
                 modelname  = modelname,
                 isPCP      = isPCP,
                 po         = po,
                 lambda     = lambda,
                 mu         = mu,
                 par        = optpar,
                 modelpar   = modelpar,
                 covmodel   = clargs,
                 Fit        = Fit)
  return(result)
}

kppmPalmLik <- function(X, Xname, po, clusters, control, weightfun, rmax, ...) {
  W <- as.owin(X)
  if(is.null(rmax))
    rmax <- rmax.rule("K", W, intensity(X))
  # identify pairs of points that contribute
  cl <- closepairs(X, rmax)
  I <- cl$i
  J <- cl$j
  dIJ <- cl$d
  # compute weights for pairs of points
  if(is.function(weightfun)) {
    wIJ <- weightfun(dIJ)
    sumweight <- sum(wIJ)
  } else {
    npairs <- length(dIJ)
    wIJ <- rep.int(1, npairs)
    sumweight <- npairs
  }
  # convert window to mask, saving other arguments for later
  dcm <- do.call.matched("as.mask",
                         append(list(w=W), list(...)),
                         sieve=TRUE)
  M         <- dcm$result
  otherargs <- dcm$otherargs
  # compute intensity at data points
  # and c.d.f. of interpoint distance in window
  if(stationary <- is.stationary(po)) {
    # stationary unmarked Poisson process
    lambda <- intensity(X)
    lambdaJ <- rep(lambda, length(J))
    # compute cdf of distance between a uniform random point in W
    # and a randomly-selected point in X 
    g <- distcdf(X, M)
    # scaling constant is (integral of intensity) * (number of points)
    gscale <- npoints(X)^2
  } else {
    # compute fitted intensity at data points and in window
    lambdaX <- fitted(po, dataonly=TRUE)
    lambda <- lambdaM <- predict(po, locations=M)
    lambdaJ <- lambdaX[J] 
    # compute cdf of distance between a uniform random point in X 
    # and a random point in W with density proportional to intensity function
    g <- distcdf(X, M, dV=lambdaM)
    # scaling constant is (integral of intensity) * (number of points)
    gscale <- integral.im(lambdaM) * npoints(X)
  }
  # trim 'g' to [0, rmax] 
  g <- g[with(g, .x) <= rmax,]
  # get pair correlation function (etc) for model
  info <- spatstatClusterModelInfo(clusters)
  pcfun      <- info$pcf
  funaux     <- info$funaux
  selfstart  <- info$selfstart
  isPCP      <- info$isPCP
  parhandler <- info$parhandler
  modelname  <- info$modelname
  # Assemble information required for computing pair correlation
  pcfunargs <- list(funaux=funaux)
  if(is.function(parhandler)) {
    # Additional parameters of cluster model are required.
    # These may be given as individual arguments,
    # or in a list called 'covmodel'
    clustargs <- if("covmodel" %in% names(otherargs))
                 otherargs[["covmodel"]] else otherargs
    clargs <- do.call(parhandler, clustargs)
    pcfunargs <- append(clargs, pcfunargs)
  } else clargs <- NULL
  # determine starting parameter values
  startpar <- selfstart(X)
  # create local function to evaluate pair correlation
  #  (with additional parameters 'pcfunargs' in its environment)
  paco <- function(d, par) {
    do.call(pcfun, append(list(par=par, rvals=d), pcfunargs))
  }
  # define objective function 
  if(!is.function(weightfun)) {
    # pack up necessary information
    objargs <- list(dIJ=dIJ, g=g, gscale=gscale,
                    sumloglam=sum(log(lambdaJ)),
                    envir=environment(paco))
    # define objective function (with 'paco' in its environment)
    # This is the log Palm likelihood
    obj <- function(par, objargs) {
      with(objargs,
           sumloglam + sum(log(paco(dIJ, par)))
           - gscale * unlist(stieltjes(paco, g, par=par)),
           enclos=objargs$envir)
    }
  } else {
    # create local function to evaluate  pair correlation(d) * weight(d)
    #  (with additional parameters 'pcfunargs', 'weightfun' in its environment)
    force(weightfun)
    wpaco <- function(d, par) {
      y <- do.call(pcfun, append(list(par=par, rvals=d), pcfunargs))
      w <- weightfun(d)
      return(y * w)
    }
    # pack up necessary information
    objargs <- list(dIJ=dIJ, wIJ=wIJ, g=g, gscale=gscale,
                    wsumloglam=sum(wIJ * log(lambdaJ)),
                    envir=environment(wpaco))
    # define objective function (with 'paco', 'wpaco' in its environment)
    # This is the log Palm likelihood
    obj <- function(par, objargs) {
      with(objargs,
           wsumloglam + sum(wIJ * log(paco(dIJ, par)))
           - gscale * unlist(stieltjes(wpaco, g, par=par)),
           enclos=objargs$envir)
    }
  }    
  # optimize it
  ctrl <- resolve.defaults(list(fnscale=-1), control, list(trace=0))
  opt <- optim(startpar, obj, objargs=objargs, control=ctrl)
  # raise warning/error if something went wrong
  signalStatus(optimStatus(opt), errors.only=TRUE)
  # meaningful model parameters
  optpar <- opt$par
  modelpar <- info$interpret(optpar, lambda)
  # infer parameter 'mu'
  if(isPCP) {
    # Poisson cluster process: extract parent intensity kappa
    kappa <- optpar[["kappa"]]
    # mu = mean cluster size
    mu <- if(stationary) lambda/kappa else eval.im(lambda/kappa)
  } else {
    # LGCP: extract variance parameter sigma2
    sigma2 <- optpar[["sigma2"]]
    # mu = mean of log intensity 
    mu <- if(stationary) log(lambda) - sigma2/2 else
          eval.im(log(lambda) - sigma2/2)    
  }
  # all info that depends on the fitting method:
  Fit <- list(method    = "palm",
              clfit     = opt,
              weightfun = weightfun,
              rmax      = rmax)
  # pack up
  result <- list(Xname      = Xname,
                 X          = X,
                 stationary = stationary,
                 clusters   = clusters,
                 modelname  = modelname,
                 isPCP      = isPCP,
                 po         = po,
                 lambda     = lambda,
                 mu         = mu,
                 par        = optpar,
                 modelpar   = modelpar,
                 covmodel   = clargs,
                 Fit        = Fit)
  return(result)
}

improve.kppm <- local({

  fnc <- function(r, eps, g){ (g(r) - 1)/(g(0) - 1) - eps}

  improve.kppm <- function(object, type=c("quasi", "wclik1", "clik1"),
                           rmax = NULL, eps.rmax = 0.01,
                           dimyx = 50, maxIter = 100, tolerance = 1e-06,
                           fast = TRUE, vcov = FALSE, fast.vcov = FALSE,
                           verbose = FALSE,
                           save.internals = FALSE) {
    verifyclass(object, "kppm")
    type <- match.arg(type)
    if(((type == "quasi" && fast) || (vcov && fast.vcov)) && !require(Matrix))
      stop(paste("Package Matrix must be installed in order for",
                 "the fast option of quasi-likelihood estimation",
                 "for cluster models to work."),
           call.=FALSE)
    gfun <- pcfmodel(object)
    X <- object$X
    win <- as.owin(X)
    ## simple (rectangular) grid quadrature scheme
    ## (using pixels with centers inside owin only)
    mask <- as.mask(win, dimyx = dimyx)
    wt <- pixellate(win, W = mask)
    wt <- wt[mask]
    Uxy <- rasterxy.mask(mask)
    U <- ppp(Uxy$x, Uxy$y, window = win, check=FALSE)
    nU <- npoints(U)
    Yu <- pixellate(X, W = mask)
    Yu <- Yu[mask]
    
    ## covariates at quadrature points
    po <- object$po
    Z <- model.images(po, mask)
    Z <- sapply(Z, "[", i=U)

    ##obtain initial beta estimate using composite likelihood
    beta0 <- coef(po)
    
    ## determining the dependence range
    if (type != "clik1" && is.null(rmax))
      {
        diamwin <- diameter(win)
        rmax <- if(fnc(diamwin, eps.rmax, gfun) >= 0) diamwin else
                uniroot(fnc, lower = 0, upper = diameter(win),
                        eps=eps.rmax, g=gfun)$root
        if(verbose)
          splat(paste0("type: ", type, ", ",
                       "dependence range: ", rmax, ", ",
                       "dimyx: ", dimyx, ", g(0) - 1:", gfun(0) -1))
      }
    ## preparing the WCL case
    if (type == "wclik1")
      Kmax <- 2*pi * integrate(function(r){r * (gfun(r) - 1)},
                               lower=0, upper=rmax)$value * exp(c(Z %*% beta0))
    ## the g()-1 matrix without tapering
    if (!fast){
      if (verbose)
        cat("computing the g(u_i,u_j)-1 matrix ...")
      gminus1 <- matrix(gfun(c(pairdist(U))) - 1, U$n, U$n)
      if (verbose)
        cat("..Done.\n")
    }
    
    if ( (fast && type == "quasi") | fast.vcov ){
      if (verbose)
        cat("computing the sparse G-1 matrix ...\n")
      ## Non-zero gminus1 entries (when using tapering)
      cp <- crosspairs(U,U,rmax)
      if (verbose)
        cat("crosspairs done\n")
      Gtap <- (gfun(cp$d) - 1)
      if(vcov){
        if(fast.vcov){
          gminus1 <- Matrix::sparseMatrix(i=cp$i, j=cp$j,
                                          x=Gtap, dims=c(U$n, U$n))
        } else{
          if(fast)
            gminus1 <- matrix(gfun(c(pairdist(U))) - 1, U$n, U$n)
        }
      }
      if (verbose & type!="quasi")
        cat("..Done.\n")
    }
       
    if (type == "quasi" && fast){
      mu0 <- exp(c(Z %*% beta0)) * wt
      mu0root <- sqrt(mu0)
      sparseG <- Matrix::sparseMatrix(i=cp$i, j=cp$j,
                                      x=mu0root[cp$i] * mu0root[cp$j] * Gtap,
                                      dims=c(U$n, U$n))
      Rroot <- Matrix::Cholesky(sparseG, perm = TRUE, Imult = 1)
      ##Imult=1 means that we add 1*I
      if (verbose)
        cat("..Done.\n")
    }
    
    ## iterative weighted least squares/Fisher scoring
    bt <- beta0
    noItr <- 1
    repeat {
      mu <- exp(c(Z %*% bt)) * wt
      mu.root <- sqrt(mu)
      ## the core of estimating equation: ff=phi
      ## in case of ql, \phi=V^{-1}D=V_\mu^{-1/2}x where (G+I)x=V_\mu^{1/2} Z
      ff <- switch(type,
                   clik1 = Z,
                   wclik1= Z/(1 + Kmax),
                   quasi = if(fast){
                     Matrix::solve(Rroot, mu.root * Z)/mu.root
                   } else{
                     solve(diag(U$n) + t(gminus1 * mu), Z)
                   }
                   )
      ##alternative
      ##R=chol(sparseG+sparseMatrix(i=c(1:U$n),j=c(1:U$n),
      ##                            x=rep(1,U$n),dims=c(U$n,U$n)))
      ##ff2 <- switch(type,
      ##              clik1 = Z,
      ##              wclik1= Z/(1 + Kmax),
      ##              quasi = if (fast)
      ##                         solve(R,solve(t(R), mu.root * Z))/mu.root
      ##                      else solve(diag(U$n) + t(gminus1 * mu), Z))
      ## print(summary(as.numeric(ff)-as.numeric(ff2)))
      ## the estimating equation: u_f(\beta)
      uf <- (Yu - mu) %*% ff
      ## inverse of minus expectation of Jacobian matrix: I_f
      Jinv <- solve(t(Z * mu) %*% ff)
      if(maxIter==0){
        ## This is a built-in early exit for vcov internal calculations
        break
      }
      deltabt <- as.numeric(uf %*% Jinv)
      if (any(!is.finite(deltabt))) {
        warning(paste("Infinite value, NA or NaN appeared",
                      "in the iterative weighted least squares algorithm.",
                      "Returning the initial intensity estimate unchanged."),
                call.=FALSE)
        return(object)
      }
      ## updating the present estimate of \beta
      bt <- bt + deltabt
      if (verbose)
        splat(paste0("itr: ", noItr, ",\nu_f: ", as.numeric(uf),
                     "\nbeta:", bt, "\ndeltabeta:", deltabt))
      if (max(abs(deltabt/bt)) <= tolerance || max(abs(uf)) <= tolerance)
        break
      if (noItr > maxIter)
        stop("Maximum number of iterations reached without convergence.")
      noItr <- noItr + 1
    }
    out <- object
    out$po$coef.orig <- beta0
    out$po$coef <- bt
    out$lambda <- predict(out$po, locations = as.mask(out$lambda))
    out$improve <- list(type = type,
                        rmax = rmax,
                        dimyx = dimyx,
                        fast = fast,
                        fast.vcov = fast.vcov)
    if(save.internals){
      out$improve <- append(out$improve, list(ff=ff, uf=uf, J.inv=Jinv))
    }
    if(vcov){
      if (verbose)
        cat("computing the asymptotic variance ...\n")
      ## variance of the estimation equation: Sigma_f = Var(u_f(bt))
      trans <- if(fast) Matrix::t else t
      Sig <- trans(ff) %*% (ff * mu) + trans(ff * mu) %*% gminus1 %*% (ff * mu)
      ## note Abdollah's G does not have mu.root inside...
      ## the asymptotic variance of \beta:
      ##         inverse of the Godambe information matrix
      out$vcov <- as.matrix(Jinv %*% Sig %*% Jinv)
    }
    return(out)
  }
  improve.kppm
})


is.kppm <- function(x) { inherits(x, "kppm")}

print.kppm <- function(x, ...) {

  isPCP <- x$isPCP
  # handle outdated objects - which were all cluster processes
  if(is.null(isPCP)) isPCP <- TRUE
  
  cat(paste(if(x$stationary) "Stationary" else "Inhomogeneous",
            if(isPCP) "cluster" else "Cox",
            "point process model\n"))

  if(nchar(x$Xname) < 20)
    cat(paste("Fitted to point pattern dataset", sQuote(x$Xname), "\n"))

  switch(x$Fit$method,
         mincon = {
           cat("Fitted by minimum contrast\n")
           cat(paste("\tSummary statistic:", x$Fit$StatName, "\n"))
         },
         clik  =,
         clik2 = {
           cat("Fitted by maximum second order composite likelihood\n")
           cat(paste("\trmax =", x$Fit$rmax, "\n"))
           if(!is.null(wtf <- x$Fit$weightfun)) {
             cat("\tweight function: ")
             print(wtf)
           }
         },
         palm = {
           cat("Fitted by maximum Palm likelihood\n")
           cat(paste("\trmax =", x$Fit$rmax, "\n"))
           if(!is.null(wtf <- x$Fit$weightfun)) {
             cat("\tweight function: ")
             print(wtf)
           }
         },
         warning(paste("Unrecognised fitting method", sQuote(x$Fit$method)))
         )

  # ............... trend .........................

  print(x$po, what="trend")

  # ..................... clusters ................
  
  cat(paste(if(isPCP) "Cluster" else "Cox",
            "model:", x$modelname, "\n"))
  cm <- x$covmodel
  if(!is.null(cm)) {
    # Covariance model - LGCP only
    cat(paste("\tCovariance model:", cm$model, "\n"))
    margs <- cm$margs
    if(!is.null(margs)) {
      nama <- names(margs)
      tags <- ifelse(nzchar(nama), paste(nama, "="), "")
      tagvalue <- paste(tags, margs)
      cat(paste("\tCovariance parameters:",
                paste(tagvalue, collapse=", "),
                "\n"))
    }
  }
  pa <- x$par
  if (!is.null(pa)) {
    cat(paste("Fitted",
              if(isPCP) "cluster" else "correlation",
              "parameters:\n"))
    print(pa)
  }

  if(!is.null(mu <- x$mu)) {
    if(isPCP) {
      cat("Mean cluster size: ")
      if(!is.im(mu)) cat(mu, "points\n") else cat("[pixel image]\n")
    } else {
      cat("Fitted mean of log of random intensity: ")
      if(!is.im(mu)) cat(mu, "\n") else cat("[pixel image]\n")
    }
  }
  invisible(NULL)
}

plot.kppm <- function(x, ..., what=c("intensity", "statistic")) {
  objectname <- short.deparse(substitute(x))
  plotem <- function(x, ..., main=dmain, dmain) { plot(x, ..., main=main) }
  what <- pickoption("plot type", what,
                    c(statistic="statistic",
                      intensity="intensity"),
                    multi=TRUE)
  # handle older objects
  Fit <- x$Fit
  if(is.null(Fit)) {
    warning("kppm object is in outdated format")
    Fit <- x
    Fit$method <- "mincon"
  } 
  inappropriate <- ((what == "intensity") & (x$stationary)) |
                   ((what == "statistic") & (Fit$method != "mincon"))
  if(any(inappropriate)) {
    what <- what[!inappropriate]
    if(length(what) == 0) return(invisible(NULL))
  }
  pauseit <- interactive() && (length(what) > 1)
  if(pauseit) opa <- par(ask=TRUE)
  for(style in what)
    switch(style,
           intensity={
             plotem(x$po, ...,
                    dmain=c(objectname, "Intensity"),
                    how="image", se=FALSE)
           },
           statistic={
             plotem(Fit$mcfit, ...,
                    dmain=c(objectname, Fit$StatName))
           })
  if(pauseit) par(opa)
  return(invisible(NULL))
}

predict.kppm <- function(object, ...) {
  predict(object$po, ...)
}

fitted.kppm <- function(object, ...) {
  fitted(object$po, ...)
}


simulate.kppm <- function(object, nsim=1, seed=NULL, ...,
                          window=NULL, covariates=NULL,
                          verbose=TRUE, retry=10) {
  starttime <- proc.time()
  verbose <- verbose && (nsim > 1)
  check.1.real(retry)
  # .... copied from simulate.lm ....
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  
  # ..................................
  # determine window for simulation results
  if(!is.null(window)) {
    stopifnot(is.owin(window))
    win <- window
  } else {
    win <- as.owin(object)
  }
  # ..................................
  # determine parameters
  mp <- as.list(object$modelpar)

  # parameter 'mu'
  # = parent intensity of cluster process
  # = mean log intensity of log-Gaussian Cox process
  
  if(is.null(covariates) && (object$stationary || is.null(window))) {
    # use existing 'mu' (scalar or image)
    mu <- object$mu
  } else {
    # recompute 'mu' using new data
    switch(object$clusters,
           Cauchy=,
           VarGamma=,
           Thomas=,
           MatClust={
             # Poisson cluster process
             kappa <- mp$kappa
             lambda <- predict(object, window=win, covariates=covariates)
             mu <- eval.im(lambda/kappa)
           },
           LGCP={
             # log-Gaussian Cox process
             sigma2 <- mp$sigma2
             lambda <- predict(object, window=win, covariates=covariates)
             mu <- eval.im(log(lambda) - sigma2/2)
           },
           stop(paste("Simulation of", sQuote(object$clusters),
                      "processes is not yet implemented"))
           )
  }

  # prepare data for execution
  out <- list()
  switch(object$clusters,
         Thomas={
           kappa <- mp$kappa
           sigma <- mp$sigma
           cmd <- expression(rThomas(kappa,sigma,mu,win))
         },
         MatClust={
           kappa <- mp$kappa
           r     <- mp$R
           cmd   <- expression(rMatClust(kappa,r,mu,win))
         },
         Cauchy = {
           kappa <- mp$kappa
           omega <- mp$omega
           cmd   <- expression(rCauchy(kappa, omega, mu, win))
         },
         VarGamma = {
           kappa  <- mp$kappa
           omega  <- mp$omega
           nu.ker <- object$covmodel$margs$nu.ker
           cmd    <- expression(rVarGamma(kappa, nu.ker, omega, mu, win))
         },
         LGCP={
           sigma2 <- mp$sigma2
           alpha  <- mp$alpha
           cm <- object$covmodel
           model <- cm$model
           margs <- cm$margs
           param <- c(0, sigma2, 0, alpha, unlist(margs))
           if(!is.im(mu)) {
             # simulate in 'win'
             cmd <- expression(rLGCP(model=model, mu=mu, param=param,
                               ..., win=win))
           } else {
             # simulate in as.owin(mu), then change window
             cmd <- expression(rLGCP(model=model, mu=mu, param=param,
                               ...)[win])
           }
         })

  # run
  if(verbose)
    cat(paste("Generating", nsim, "simulations... "))
  for(i in 1:nsim) {
    out[[i]] <- try(eval(cmd))
    if(verbose) progressreport(i, nsim)
  }
  # detect failures
  if(any(bad <- unlist(lapply(out, inherits, what="try-error")))) {
    nbad <- sum(bad)
    gripe <- paste(nbad,
                   ngettext(nbad, "simulation was", "simulations were"),
                   "unsuccessful")
    if(verbose) cat(paste(gripe, "\n"))
    if(retry <= 0) {
      fate <- "returned as NULL"
      out[bad] <- list(NULL)
    } else {
      if(verbose) cat("Retrying...")
      ntried <- 0
      while(ntried < retry) {
        ntried <- ntried + 1
        for(j in which(bad))
          out[[j]] <- try(eval(cmd))
        bad <- unlist(lapply(out, inherits, what="try-error"))
        nbad <- sum(bad)
        if(nbad == 0) break
      }
      if(verbose) cat("Done.\n")
      fate <- if(nbad == 0) "all recomputed" else
              paste(nbad, "simulations still unsuccessful")
      fate <- paste(fate, "after", ntried,
                    ngettext(ntried, "further try", "further tries"))
    }
    warning(paste(gripe, fate, sep=": "))
  }
  if(verbose)
    cat("Done.\n")
  # pack up
  out <- as.listof(out)
  names(out) <- paste("Simulation", 1:nsim)
  attr(out, "seed") <- RNGstate
  out <- timed(out, starttime=starttime)
  return(out)
}

formula.kppm <- function(x, ...) {
  formula(x$po, ...)
}

terms.kppm <- function(x, ...) {
  terms(x$po, ...)
}

labels.kppm <- function(object, ...) {
  labels(object$po, ...)
}

update.kppm <- function(object, trend=~1, ..., clusters=NULL) {
  if(!missing(trend))
    trend <- update(formula(object), trend)
  if(is.null(clusters))
    clusters <- object$clusters
  out <- do.call(kppm,
                 resolve.defaults(list(trend=trend, clusters=clusters),
                                  list(...),
                                  list(X=object$X)))
  out$Xname <- object$Xname
  return(out)
}

unitname.kppm <- function(x) {
  return(unitname(x$X))
}

"unitname<-.kppm" <- function(x, value) {
  unitname(x$X) <- value
  if(!is.null(x$Fit$mcfit)) {
    unitname(x$Fit$mcfit) <- value
  } else if(is.null(x$Fit)) {
    warning("kppm object in outdated format")
    if(!is.null(x$mcfit))
      unitname(x$mcfit) <- value
  }
  return(x)
}

as.fv.kppm <- function(x) as.fv(x$Fit$mcfit)

coef.kppm <- function(object, ...) {
  return(coef(object$po))
}


Kmodel.kppm <- function(model, ...) {
  Kpcf.kppm(model, what="K")
}

pcfmodel.kppm <- function(model, ...) {
  Kpcf.kppm(model, what="pcf")
}

Kpcf.kppm <- function(model, what=c("K", "pcf")) {
  what <- match.arg(what)
  # Extract function definition from internal table
  clusters <- model$clusters
  tableentry <- spatstatClusterModelInfo(clusters)
  if(is.null(tableentry))
    stop("No information available for", sQuote(clusters), "cluster model")
  fun <- tableentry[[what]]
  if(is.null(fun))
    stop("No expression available for", what, "for", sQuote(clusters),
         "cluster model")
  # Extract model parameters
  par <- model$par
  # Extract auxiliary definitions (if applicable)
  funaux <- tableentry$funaux
  # Extract covariance model (if applicable)
  cm <- model$covmodel
  model <- cm$model
  margs <- cm$margs
  #
  f <- function(r) as.numeric(fun(par=par, rvals=r,
                                  funaux=funaux, model=model, margs=margs))
  return(f)
}

is.stationary.kppm <- function(x) {
  return(x$stationary)
}

is.poisson.kppm <- function(x) {
  switch(x$clusters,
         Cauchy=,
         VarGamma=,
         Thomas=,
         MatClust={
           # Poisson cluster process
           mu <- x$mu
           return(!is.null(mu) && (max(mu) == 0))
         },
         LGCP = {
           # log-Gaussian Cox process
           sigma2 <- x$par[["sigma2"]]
           return(sigma2 == 0)
         },
         return(FALSE))
}

# extract ppm component

as.ppm.kppm <- function(object) {
  object$po
}

# other methods that pass through to 'ppm'

as.owin.kppm <- function(W, ..., from=c("points", "covariates"), fatal=TRUE) {
  from <- match.arg(from)
  as.owin(as.ppm(W), ..., from=from, fatal=fatal)
}

domain.kppm <- Window.kppm <- function(X, ..., from=c("points", "covariates")) {
  from <- match.arg(from)
  as.owin(X, from=from)
}

model.images.kppm <- function(object, W=as.owin(object), ...) {
  model.images(as.ppm(object), W=W, ...)
}

model.matrix.kppm <- function(object, data=model.frame(object), ...,
                              keepNA=TRUE) {
  model.matrix(as.ppm(object), data=data, ..., keepNA=keepNA)
}

model.frame.kppm <- function(formula, ...) {
  model.frame(as.ppm(formula), ...)
}

