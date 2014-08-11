#
#	$Revision: 1.32 $	$Date: 2013/07/19 03:55:36 $
#
#    ppm()
#          Fit a point process model to a two-dimensional point pattern
#
#

"ppm" <- 
function(Q,
         trend = ~1,
	 interaction = Poisson(),
         ..., 
         covariates = NULL,
         covfunargs = list(),
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
  
  if(!(method %in% c("mpl", "ho", "logi")))
    stop(paste("Unrecognised fitting method", sQuote(method)))
  cl <- match.call()
  if(is.null(callstring)) 
    callstring <- paste(short.deparse(sys.call()), collapse="")

  if(is.ppp(Q) && is.marked(Q) && !is.multitype(Q)) 
    stop(paste("ppm is not yet implemented for marked point patterns,",
               "other than multitype patterns."))
  if(!(is.ppp(Q) || inherits(Q, "quad")))
    stop("Argument Q must be a point pattern or a quadrature scheme")

# Ensure interaction is fully defined  
  if(is.null(interaction)) 
    interaction <- Poisson()
  if(!is.null(ss <- interaction$selfstart)) {
    # invoke selfstart mechanism to fix all parameters
    X <- if(is.ppp(Q)) Q else Q$data
    interaction <- ss(X, interaction)
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
    too.large <- infin || (eroded.areas(as.owin(Q), rbord) == 0)
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

