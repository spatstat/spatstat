#
#  logistic.R
#
#   $Revision: 1.19 $  $Date: 2015/04/02 02:17:19 $
#
#  Logistic likelihood method - under development
#

logi.engine <- function(Q,
                        trend = ~1,
                        interaction,
                        ...,
                        covariates=NULL,
                        subsetexpr=NULL,
                        correction="border",
                        rbord=reach(interaction),
                        covfunargs=list(),
                        allcovar=FALSE,
                        vnamebase=c("Interaction", "Interact."),
                        vnameprefix=NULL,
                        justQ = FALSE,
                        savecomputed = FALSE,
                        precomputed = NULL,
                        VB=FALSE
                        ){
  if(is.null(trend)) trend <- ~1 
  if(is.null(interaction)) interaction <- Poisson()
  want.trend <- !identical.formulae(trend, ~1)
  want.inter <- !is.poisson(interaction)
  want.subset <- !is.null(subsetexpr)
  # validate choice of edge correction
  correction <- pickoption("correction", correction,
                           c(border="border",
                             periodic="periodic",
                             isotropic="isotropic",
                             Ripley="isotropic",
                             trans="translate",
                             translate="translate",
                             translation="translate",
                             none="none"))
  # rbord applies only to border correction
  if(correction == "border") {
    check.1.real(rbord, "In ppm")
    explain.ifnot(rbord >= 0, "In ppm")
  } else rbord <- 0
  # backdoor stuff
  if(!missing(vnamebase)) {
    if(length(vnamebase) == 1)
      vnamebase <- rep.int(vnamebase, 2)
    if(!is.character(vnamebase) || length(vnamebase) != 2)
      stop("Internal error: illegal format of vnamebase")
  }
  if(!is.null(vnameprefix)) {
    if(!is.character(vnameprefix) || length(vnameprefix) != 1)
      stop("Internal error: illegal format of vnameprefix")
  }
  # create dummy points
  if(inherits(Q, "ppp")){
    Xplus <- Q
    Q <- quadscheme.logi(Xplus, ...)
    D <- Q$dummy
    Dinfo <- Q$param
  } else if(checkfields(Q, c("data", "dummy"))) {
    Xplus <- Q$data
    D <- Q$dummy
    Dinfo <- Q$param
    if(is.null(Dinfo)){
      Dinfo <- list(how="given", rho=npoints(D)/(area(D)*markspace.integral(D)))
    }
    Q <- quadscheme.logi(Xplus, D)
  } else stop("Format of object Q is not understood")
  if (justQ) 
    return(Q)
  ### Dirty way of recording arguments so that the model can be refitted later (should probably be done using call, eval, envir, etc.):
  extraargs <- list(covfunargs = covfunargs, allcovar = allcovar, vnamebase = vnamebase, vnameprefix = vnameprefix)
  extraargs <- append(extraargs, list(...))
  ## Dummy intensity
  if(correction == "border" && Dinfo$how=="grid"){
    Dbord <- D[bdist.points(D)>=rbord]
    Dinfo$rho <- npoints(Dbord)/(eroded.areas(as.owin(Dbord), rbord)*markspace.integral(Dbord))
  }
  rho <- Dinfo$rho
  ##Setting the B from Barker dynamics (relative to dummy intensity)
  B <- list(...)$Barker
  if(is.null(B))
    B <- 1
  B <- B*rho
  Dinfo <- append(Dinfo, list(B=B))
  Dinfo <- append(Dinfo, list(extraargs=extraargs))
  # 
  Wplus <- as.owin(Xplus)
  nXplus <- npoints(Xplus)
  U <- superimpose(Xplus, D, W=Wplus, check=FALSE)
#  E <- equalpairs(U, Xplus, marked = is.marked(Xplus))
  E <- cbind(1:nXplus, 1:nXplus)
#  
  computed <- if (savecomputed) list(X = Xplus, Q = Q, U = U) else list()
  # assemble covariate data frame
  if(want.trend || want.subset) {
    tvars <- variablesinformula(trend)
    if(want.subset)
      tvars <- union(tvars, all.vars(subsetexpr))
    ## resolve 'external' covariates
    externalvars <- setdiff(tvars, c("x", "y", "marks"))
    tenv <- environment(trend)
    covariates <- getdataobjects(externalvars, tenv, covariates, fatal=TRUE)
    ## 
    wantxy <- c("x", "y") %in% tvars
    wantxy <- wantxy | rep.int(allcovar, 2)
    cvdf <- data.frame(x=U$x, y=U$y)[, wantxy, drop=FALSE]
    if(!is.null(covariates)) {
      df <- mpl.get.covariates(covariates, U, "quadrature points", covfunargs)
      cvdf <- cbind(cvdf, df)
    }
    wantmarks <- "marks" %in% tvars
    if(wantmarks) cvdf <- cbind(cvdf, marks = marks(U))
  } else cvdf <- NULL
  # evaluate interaction sufficient statistics
  if (!is.null(ss <- interaction$selfstart)) 
    interaction <- ss(Xplus, interaction)
  V <- evalInteraction(Xplus, U, E, interaction, correction, precomputed = precomputed, savecomputed = savecomputed)
  if(!is.matrix(V))
    stop("evalInteraction did not return a matrix")
  if (savecomputed) 
    computed <- append(computed, attr(V, "computed"))
  IsOffset <- attr(V, "IsOffset")
  if(is.null(IsOffset)) IsOffset <- rep.int(FALSE, ncol(V))
  # determine names
  if(ncol(V) > 0) {
    Vnames <- colnames(V)
    if(is.null(Vnames)) {
      nc <- ncol(V)
      Vnames <- if(nc == 1) vnamebase[1] else paste(vnamebase[2], 1:nc, sep="")
      colnames(V) <- Vnames
    } else if(!is.null(vnameprefix)) {
      Vnames <- paste(vnameprefix, Vnames, sep="")
      colnames(V) <- Vnames
    }
  } else Vnames <- character(0)
  # combine all data
  glmdata <- as.data.frame(V)
  if(!is.null(cvdf)) glmdata <- cbind(glmdata, cvdf)
  # construct response and weights
  ok <- if(correction == "border") (bdist.points(U) >= rbord) else rep.int(TRUE, npoints(U))
  # Keep only those quadrature points for which the
  # conditional intensity is nonzero.
  KEEP  <- if(ncol(V)>0) matrowall(V != -Inf) else rep.int(TRUE, npoints(U))
  ok <- ok & KEEP
  wei <- c(rep.int(1,npoints(Xplus)),rep.int(B/rho,npoints(D)))
  resp <- c(rep.int(1,npoints(Xplus)),rep.int(0,npoints(D)))
  ## User-defined subset:
  if(!is.null(subsetexpr)) {
    USERSUBSET <- eval(subsetexpr, glmdata, environment(trend))
    ok <- ok & USERSUBSET
  }
  # add offset, subset and weights to data frame
  # using reserved names beginning with ".logi."
  glmdata <- cbind(glmdata,
                   .logi.Y = resp,
                   .logi.B = B,
                   .logi.w = wei,
                   .logi.ok =ok)
  # build glm formula 
  # (reserved names begin with ".logi.")
  trendpart <- paste(as.character(trend), collapse=" ")
  fmla <- paste(".logi.Y ", trendpart)
  # Interaction terms
  if(want.inter) {
    VN <- Vnames
    # enclose offset potentials in 'offset(.)'
    if(any(IsOffset))
      VN[IsOffset] <- paste("offset(", VN[IsOffset], ")", sep="")
    fmla <- paste(c(fmla, VN), collapse="+")
  }
  # add offset intrinsic to logistic technique
  fmla <- paste(fmla, "offset(-log(.logi.B))", sep="+")
  fmla <- as.formula(fmla)
  # to satisfy package checker: 
  .logi.B <- B
  .logi.w <- wei
  .logi.ok  <- ok
  .logi.Y   <- resp
  # go
  ##fit <- glm(fmla, data=glmdata,
  ##           family=binomial(), subset = .logi.ok, weights = .logi.w)
  fit <- if(VB) 
           vblogit.fmla(fmla, data = glmdata, 
                        subset = .logi.ok, weights = .logi.w, ...)
         else 
           glm(fmla, data = glmdata, 
               family = binomial(), subset = .logi.ok, weights = .logi.w)
  environment(fit$terms) <- sys.frame(sys.nframe())
  ## Fitted coeffs
  co <- coef(fit)
  fitin <- fii(interaction, co, Vnames, IsOffset)

  ## Max. value of log-likelihood:
  maxlogpl <- logLik(fit) + sum(ok*resp*log(B))

  # Stamp with spatstat version number
  spv <- package_version(versionstring.spatstat())
  the.version <- list(major=spv$major,
                      minor=spv$minor,
                      release=spv$patchlevel,
                      date="$Date: 2015/04/02 02:17:19 $")

  ## Compile results
  fit <- list(method      = "logi",
              fitter      = "glm",
              projected   = FALSE,
              coef        = co,
              trend       = trend,
              interaction = interaction,
              Q           = Q,
              correction  = correction,
              rbord       = rbord,
              terms       = terms(trend),
              version     = the.version,
              fitin       = fitin,
              maxlogpl    = maxlogpl,
              covariates  = mpl.usable(covariates),
#              varcov      = if(VB) fit$S else NULL,
              internal    = list(Vnames  = Vnames,
                                 IsOffset=IsOffset,
                                 glmdata = glmdata,
                                 glmfit = fit,
                                 logistic = Dinfo,
                                 computed = computed,
                                 VB = if(VB) TRUE else NULL,
                                 priors = if(VB) fit$priors else NULL
                                 )
              )
  class(fit) <- "ppm"
  return(fit)
}


forbid.logi <- function(object) {
  if(object$method == "logi")
    stop("Sorry, this is not implemented for method=\'logi\'")
  return(invisible(NULL))
}

logi.dummy <- function(X, dummytype = "stratrand", nd = NULL, mark.repeat = FALSE, ...){
  ## Resolving nd inspired by default.n.tiling
  if(is.null(nd)){
    nd <- spatstat.options("ndummy.min")
    if(inherits(X, "ppp"))
      nd <- pmax(nd, 10 * ceiling(2 * sqrt(X$n)/10))
  }
  nd <- ensure2vector(nd)
  marx <- is.multitype(X)
  if(marx)
    lev <- levels(marks(X))
  if(marx && mark.repeat){
    N <- length(lev)
    Dlist <- inDlist <- vector("list", N)
  } else{
    N <- 1
  }
  W <- as.owin(X)
  type <- match.arg(dummytype, c("stratrand", "binomial", "poisson", "grid", "transgrid"))
  B <- boundingbox(W)
  rho <- nd[1]*nd[2]/area(B)
  Dinfo <- list(nd=nd, rho=rho, how=type)
  ## Repeating dummy process for each mark type 1:N (only once if unmarked or mark.repeat = FALSE)
  for(i in 1:N){
    switch(type,
           stratrand={
             D <- as.ppp(stratrand(B, nd[1], nd[2]), W = B)
             inD <- which(inside.owin(D, w = W))
             D <- D[W]
             inD <- paste(i,inD,sep="_")
           },
           binomial={
             D <- runifpoint(nd[1]*nd[2], win=B)
             D <- D[W]
           },
           poisson={
             D <- rpoispp(rho, win = W)
           },
           grid={
             D <- as.ppp(gridcenters(B, nd[1], nd[2]), W = B)
             inD <- which(inside.owin(D, w = W))
             D <- D[W]
             inD <- paste(i,inD,sep="_")
           },
           transgrid={
             D <- as.ppp(gridcenters(B, nd[1], nd[2]), W = B)
             dxy <- c(diff(D$window$xrange),diff(D$window$yrange))/(2*nd)
             coords(D) <- coords(D)+matrix(runif(2,-dxy,dxy),npoints(D),2,byrow=TRUE)
             inD <- which(inside.owin(D, w = W))
             D <- D[W]
             inD <- paste(i,inD,sep="_")
           },
         stop("unknown dummy type"))
    if(marx && mark.repeat){
      marks(D) <- factor(lev[i], levels = lev)
      Dlist[[i]] <- D
      if(type %in% c("stratrand","grid","transgrid"))
        inDlist[[i]] <- inD
    }
  }
  if(marx && mark.repeat){
    inD <- Reduce(append, inDlist)
    D <- Reduce(superimpose, Dlist)
  }
  if(type %in% c("stratrand","grid","transgrid"))
    Dinfo <- append(Dinfo, list(inD=inD))
  if(marx && !mark.repeat){
    marks(D) <- sample(factor(lev, levels=lev), npoints(D), replace = TRUE)
    Dinfo$rho <- Dinfo$rho/length(lev)
  }
  attr(D, "dummy.parameters") <- Dinfo
  return(D)
}

quadscheme.logi <- function(data, dummy, dummytype = "stratrand", nd = NULL, mark.repeat = FALSE, ...){
  data <- as.ppp(data)
  ## If dummy is missing we generate dummy pattern with logi.dummy.
  if(missing(dummy))
    dummy <- logi.dummy(data, dummytype, nd, mark.repeat, ...)
  Dinfo <- attr(dummy, "dummy.parameters")
  D <- as.ppp(dummy)
  if(is.null(Dinfo))
    Dinfo <- list(how="given", rho=npoints(D)/(area(D)*markspace.integral(D)))
  ## Weights:
  n <- npoints(data)+npoints(D)
  w <- area(Window(data))/n
  Q <- quad(data, D, rep(w,n), param=Dinfo)
  class(Q) <- c("logiquad", class(Q))
  return(Q)
}

summary.logiquad <- function(object, ..., checkdup=FALSE) {
  verifyclass(object, "logiquad")
  s <- list(
       data  = summary.ppp(object$data, checkdup=checkdup),
       dummy = summary.ppp(object$dummy, checkdup=checkdup),
       param = object$param)
  class(s) <- "summary.logiquad"
  return(s)
}

print.summary.logiquad <- function(x, ..., dp=3) {
  cat("Quadrature scheme = data + dummy\n")
  Dinfo <- x$param
  if(is.null(Dinfo))
    cat("created by an unknown function.\n")
  cat("Data pattern:\n")
  print(x$data, dp=dp)

  cat("\n\nDummy pattern:\n")
  # How they were computed
    switch(Dinfo$how,
           stratrand={
             cat(paste("(Stratified random dummy points,",
                       paste(Dinfo$nd, collapse=" x "),
                       "grid of cells)\n"))
           },
           binomial={
             cat("(Binomial dummy points)\n")
           },
           poisson={
             cat("(Poisson dummy points)\n")
           },
           grid={
             cat(paste("(Fixed grid of dummy points,",
                       paste(Dinfo$nd, collapse=" x "),
                       "grid)\n"))
           },
           transgrid={
             cat(paste("(Random translation of fixed grid of dummy points,",
                       paste(Dinfo$nd, collapse=" x "),
                       "grid)\n"))
           },
           given=cat("(Dummy points given by user)\n")
       )
  # Description of them
  print(x$dummy, dp=dp)

  return(invisible(NULL))
}
