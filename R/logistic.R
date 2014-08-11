#
#  logistic.R
#
#   $Revision: 1.6 $  $Date: 2013/04/25 06:37:43 $
#
#  Logistic likelihood method - under development
#

logi.engine <- local({
  # The true syntax is hidden 
  visible.logi.engine <- function(...) actual.logi.engine(...)
  
  actual.logi.engine <- function(Q, trend = ~1, interaction, ...,
                                 covariates=NULL,
                                 nd = 10,
                                 correction="border", rbord=reach(interaction),
                                 covfunargs=list(),
                                 allcovar=FALSE,
                                 vnamebase=c("Interaction", "Interaction."),
                                 vnameprefix=NULL
                        ){
  if(is.null(trend)) trend <- ~1 
  if(is.null(interaction)) interaction <- Poisson()
  want.trend <- !identical.formulae(trend, ~1)
  want.inter <- !is.poisson(interaction)
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
    nd <- ensure2vector(nd)
    D <- as.ppp(stratrand(as.owin(Xplus), nd[1], nd[2]), W = as.owin(Xplus))
    Q <- quad(Xplus, D)
    Dinfo <- list(how="stratrand", args=nd)
  } else if(checkfields(Q, c("data", "dummy"))) {
    Xplus <- Q$data
    D <- Q$dummy
    Dinfo <- list(how="given")
  } else stop("Format of object Q is not understood")
  rho <- npoints(D)/area.owin(D)
  B <- rho ##Setting the B from Barker dynamics
  Dinfo <- append(Dinfo, list(B=B))
  # 
  Wplus <- as.owin(Xplus)
  W <- erosion(Wplus,rbord)
  D <- D[W]
  U <- superimpose(Xplus, D)
  E <- equalpairs(U, Xplus)
  # assemble covariate data frame
  if(want.trend) {
    tvars <- variablesinformula(trend)
    wantxy <- c("x", "y") %in% tvars
    wantxy <- wantxy | rep.int(allcovar, 2)
    cvdf <- data.frame(x=U$x, y=U$y)[, wantxy, drop=FALSE]
    if(!is.null(covariates)) {
      df <- mpl.get.covariates(covariates, U, "quadrature points", covfunargs)
      cvdf <- cbind(cvdf, df)
    }
  } else cvdf <- NULL
  # evaluate interaction sufficient statistics
  V <- evalInteraction(Xplus, U, E, interaction, correction)
  if(!is.matrix(V))
    stop("evalInteraction did not return a matrix")
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
  ok <- if(correction == "border") inside.owin(U,,W) else rep.int(TRUE, npoints(U))
  wei <- c(rep.int(1,npoints(Xplus)),rep.int(B/rho,npoints(D)))
  resp <- c(rep.int(1,npoints(Xplus)),rep.int(0,npoints(D)))
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
  fit <- glm(fmla, data=glmdata,
             family=binomial(), subset = .logi.ok, weights = .logi.w)
  ## Fitted coeffs
  co <- coef(fit)
  fitin <- fii(interaction, co, Vnames, IsOffset)
  ## Compile results
  fit <- list(method      = "logi",
              fitter      = "glm",
              coef        = co,
              trend       = trend,
              interaction = interaction,
              Q           = Q,
              correction  = correction,
              rbord       = rbord,
              version     = versionstring.spatstat(),
              fitin       = fitin,
              internal    = list(Vnames  = Vnames,
                                 IsOffset=IsOffset,
                                 logistic=Dinfo)
              )
  class(fit) <- c("logippm", "ppm")
  return(fit)
}

  visible.logi.engine
})

print.logippm <- function(x, ...) {
  cat("Point process model fitted by logistic regression\n")
  if(identical.formulae(x$trend, ~1)) {
    cat("Stationary ")
    print(fitin(x))
  } else {
    cat(paste("Trend formula:", paste(x$trend, collapse=" "), "\n"))
    cat("Interaction: ")
    print(fitin(x))
  }
  cat("Fitted coefficients:\n")
  print(coef(x))
  if(!is.null(Dinfo <- x$internal$logistic)) {
    cat("\n")
    switch(Dinfo$how,
           stratrand={
             cat(paste("Stratified random dummy points,",
                       paste(Dinfo$args, collapse=" x "),
                       "grid of cells\n"))
           },
           given=cat("Dummy points given by user\n"),
           warning("Unrecognised format of internal data"))
    cat(paste("Dummy intensity", Dinfo$B, "\n"))
  }
  return(invisible(NULL))
}

summary.logippm <- function(object, ...) {
  print(object, ...)
}

forbid.logi <- function(object) {
  if(object$method == "logi")
    stop("Sorry, this is not implemented for method=\'logi\'")
  return(invisible(NULL))
}

