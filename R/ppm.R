#
#	$Revision: 1.49 $	$Date: 2014/11/10 11:27:55 $
#
#    ppm()
#          Fit a point process model to a two-dimensional point pattern
#
#

ppm <- function(Q, ...) {
  UseMethod("ppm")
}


ppm.formula <- function(Q, interaction=NULL, ..., data=NULL, subset) {
  ## remember call
  callstring <- short.deparse(sys.call())
  cl <- match.call()

  ########### INTERPRET FORMULA ##############################
  
  if(!inherits(Q, "formula"))
    stop(paste("Argument 'Q' should be a formula"))
  formula <- Q
  
  ## check formula has LHS and RHS. Extract them
  if(length(formula) < 3)
    stop(paste("Formula must have a left hand side"))
  Yexpr <- formula[[2]]
  trend <- formula[c(1,3)]
  
  ## FIT #######################################
  thecall <- if(missing(subset)) {
    call("ppm", Q=Yexpr, trend=trend, data=data, interaction=interaction)
  } else {
    call("ppm", Q=Yexpr, trend=trend, data=data, interaction=interaction,
         subset=substitute(subset))
  }
  ncall <- length(thecall)
  argh <- list(...)
  nargh <- length(argh)
  if(nargh > 0) {
    thecall[ncall + 1:nargh] <- argh
    names(thecall)[ncall + 1:nargh] <- names(argh)
  }
  callenv <- parent.frame()
  if(!is.null(data)) callenv <- list2env(data, parent=callenv)
  result <- eval(thecall, envir=callenv)

  result$call <- cl
  result$callstring <- callstring
  result$callframe <- parent.frame()
  
  return(result)
}


ppm.quad <- ppm.ppp <- ppm.default <- 
function(Q,
         trend = ~1,
	 interaction = Poisson(),
         ..., 
         covariates = data,
         data = NULL,
         covfunargs = list(),
         subset,
	 correction="border",
	 rbord = reach(interaction),
         use.gam=FALSE,
         method = "mpl",
         forcefit=FALSE,
         project=FALSE,
         nd = NULL,
         eps = NULL,
         gcontrol=list(),
         nsim=100,
         nrmh=1e5,
         start=NULL,
         control=list(nrep=nrmh),
         verb=TRUE,
         callstring=NULL
) {
  Qname <- short.deparse(substitute(Q))

  subsetexpr <- if(!missing(subset)) substitute(subset) else NULL

  if(!(method %in% c("mpl", "ho", "logi")))
    stop(paste("Unrecognised fitting method", sQuote(method)))
  if(inherits(Q, "logiquad")){
    if(missing(method))
      method <- "logi"
    if(method != "logi")
      stop(paste("Only method =", sQuote("logi"),
                 "makes sense when Q is of type", sQuote("logiquad")))
  }
  cl <- match.call()
  if(is.null(callstring)) 
    callstring <- paste(short.deparse(sys.call()), collapse="")

  if(is.ppp(Q) && is.marked(Q) && !is.multitype(Q)) 
    stop(paste("ppm is not yet implemented for marked point patterns,",
               "other than multitype patterns."))
  if(!(is.ppp(Q) ||
       inherits(Q, "quad") ||
       checkfields(Q, c("data", "dummy")))) {
    stop("Argument Q must be a point pattern or a quadrature scheme")
  }
  X <- if(is.ppp(Q)) Q else Q$data

  ## Ensure interaction is fully defined  
  if(is.null(interaction)) 
    interaction <- Poisson()
  if(!is.null(ss <- interaction$selfstart)) {
    # invoke selfstart mechanism to fix all parameters
    interaction <- ss(X, interaction)
  }

  if(inherits(trend, "formula")) {
    ## handle "." in formula, representing all variables in 'data'
    if("." %in% variablesinformula(trend)) {
      if(is.null(covariates))
        stop("Cannot expand '.' since 'data' is not present", call.=FALSE)
      rhs <- paste(names(covariates), collapse=" + ")
      allmaineffects <- as.formula(paste("~", rhs))
      environment(allmaineffects) <- environment(trend)
      trend <- update(allmaineffects, trend)
    }
    ## expand polynom() in formula
    if(spatstat.options("expand.polynom"))
      trend <- expand.polynom(trend)
  }
  
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
  
  # validate rbord 
  if(correction == "border") {
    # rbord for border correction
    rbord.given <- !missing(rbord) && !is.null(rbord)
    if(is.null(rbord))
      rbord <- reach(interaction)
    infin <- is.infinite(rbord)
    too.large <- infin || (eroded.areas(as.owin(X), rbord) == 0)
    if(too.large) {
      whinge <-
        paste(if(rbord.given) "rbord" else "the reach of this interaction",
              if(infin) "is infinite or unknown;"
              else "is too large for this window;",
              "please specify",
              if(rbord.given) "a smaller value of",
              "rbord, or use a different edge correction")
      stop(whinge)
    }
  } else {
    # rbord must be numeric to satisfy mpl.engine
    if(is.null(rbord))
      rbord <- 0
  }

  if(method == "logi") {
    fitLOGI <- logi.engine(Q=Q, trend=trend,
                           interaction=interaction,
                           covariates=covariates,
                           covfunargs=covfunargs,
                           subsetexpr=subsetexpr,
                           correction=correction,
                           rbord=rbord,
                           use.gam=use.gam,
                           forcefit=forcefit,
                           nd = nd,
                           gcontrol=gcontrol,
                           callstring=callstring,
                           ...)
    fitLOGI$Qname <- Qname
    fitLOGI$call <- cl
    fitLOGI$callstring <- callstring
    fitLOGI$callframe <- parent.frame()
    if(project && !valid.ppm(fitLOGI))
      fitLOGI <- project.ppm(fitLOGI)
    return(fitLOGI)
  }
  
  # fit by maximum pseudolikelihood
  fitMPL <- mpl.engine(Q=Q, trend=trend,
                       interaction=interaction,
                       covariates=covariates,
                       covfunargs=covfunargs,
                       subsetexpr=subsetexpr,
                       correction=correction,
                       rbord=rbord,
                       use.gam=use.gam,
                       forcefit=forcefit,
                       nd = nd,
                       eps = eps, 
                       gcontrol=gcontrol,
                       callstring=callstring,
                       ...)
  fitMPL$Qname <- Qname

  if(!is.ppm(fitMPL)) {
    # internal use only - returns some other data
    return(fitMPL)
  }
  
  fitMPL$call <- cl
  fitMPL$callstring <- callstring
  fitMPL$callframe <- parent.frame()

  if(project && !valid.ppm(fitMPL))
    fitMPL <- project.ppm(fitMPL)
  
  if(method == "mpl" || is.poisson.ppm(fitMPL))
    return(fitMPL)

  fitHO <- ho.engine(fitMPL, nsim=nsim, nrmh=nrmh, start=start,
                     control=control, verb=verb)

  if(is.null(fitHO))
    return(fitMPL)
  
  if(project && !valid.ppm(fitHO))
    fitHO <- project.ppm(fitHO)
  
  return(fitHO)
}

