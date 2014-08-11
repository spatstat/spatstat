#
# ippm.R
#
#   $Revision: 2.11 $   $Date: 2014/04/05 08:09:48 $
#
# Fisher scoring algorithm for irregular parameters in ppm trend
#

ippm <- local({
  
  ippm <- function(...,
                   iScore=NULL, 
                   start=list(),
                   covfunargs=start,
                   nlm.args=list(),
                   silent=FALSE) {
    ## validate
    if(!is.list(start) || length(start) == 0)
      stop("start should be a list of initial values for irregular parameters")
    if(!is.null(iScore)) {
      if(!is.list(iScore) || length(iScore) != length(start))
        stop("iScore should be a list of the same length as start")
      stopifnot(identical(names(iScore), names(start)))
      if(!all(unlist(lapply(iScore, is.function))))
        stop("iScore should be a list of functions")
    }
    ##
    smap <- match(names(start), names(covfunargs))
    if(any(is.na(smap)))
      stop("variables in start should be a subset of variables in covfunargs")
    covfunargs[smap] <- start
    ## fit the initial model and extract information
    fit0 <- ppm(..., covfunargs=covfunargs)
    lpl0 <- fit0$maxlogpl
    p <- length(coef(fit0))
    ## examine covariates and trend
    covariates <- fit0$covariates
    isfun <- unlist(lapply(covariates, is.function))
    covfuns <- covariates[isfun]
    ## determine which covariates depend on which irregular parameters
    pnames <- names(start)
    hasarg <- function(f,a) { a %in% names(formals(f)) }
    depmat <- matrix(FALSE, nrow=length(covfuns), ncol=length(pnames))
    rownames(depmat) <- names(covfuns)
    colnames(depmat) <- pnames
    for(j in 1:length(pnames))
      depmat[,j] <- unlist(lapply(covfuns, hasarg, pnames[j]))
    ## find covariates that depend on ANY irregular parameter 
    depvar <- rownames(depmat)[apply(depmat, 1, any)]
    ## check that these covariates appear only in offset terms
    covnames.fitted <- model.covariates(fit0, fitted=TRUE,  offset=FALSE)
    if(any(uhoh <- depvar %in% covnames.fitted))
      stop(paste(ngettext(sum(uhoh), "The covariate", "The covariates"),
                 commasep(sQuote(depvar[uhoh])),
                 "should appear only in offset terms"))
    ## check that every irregular parameter to be updated appears somewhere 
    cov.names.offset <- model.covariates(fit0, fitted=FALSE,  offset=TRUE)
    covfun.names.offset <- intersect(cov.names.offset, names(covfuns))
    usearg <- apply(depmat[covfun.names.offset, , drop=FALSE], 2, any)
    if(!all(usearg)) {
      nbad <- sum(!usearg)
      warning(paste("Cannot maximise over the irregular",
                    ngettext(nbad, "parameter", "parameters"),
                    commasep(sQuote(names(usearg)[!usearg])),
                    ngettext(nbad, "because it is", "because they are"),
                    "not used in any term of the model"))
      ## restrict 
      start <- start[usearg]
      if(!is.null(iScore)) iScore <- iScore[usearg]
      pnames <- names(start)
    }
    ## define objective function
    fdata <- list(fit0=fit0,
                  nreg=length(coef(fit0)),
                  covfunargs=covfunargs,
                  smap=smap,
                  pnames=pnames,
                  iScore=iScore)
    f <- function(param, thedata=fdata) {
      with(thedata, {
        ## fit model with current irregular parameters
        param <- as.list(param)
        names(param) <- pnames
        covfunargs[smap] <- param
        fit <- update(fit0, covfunargs=covfunargs, use.internal=TRUE)
        lpl <- logLik(fit, warn=FALSE)
        ## return negative logL because nlm performs *minimisation*
        value <- -as.numeric(lpl)
        ## compute derivatives
        stuff <- ppmInfluence(fit, what="derivatives",
                               iScore=iScore,
                               iArgs=param)
        score <- stuff$deriv$score
        if(length(score) == length(coef(fit)) + length(param)) 
          attr(value, "gradient") <- -score[-(1:nreg), drop=FALSE]
        ## attr(value, "hessian") <- -hess[-(1:nreg), -(1:nreg), drop=FALSE]
        return(value)
      })
    }
    ## minimise objective
    startvec <- unlist(start)
    typsize <- abs(startvec)
    typsize <- pmax(typsize, min(typsize[typsize > 0]))
    g <- do.call("nlm",
                 resolve.defaults(list(f=f, p=startvec, thedata=fdata),
                                  nlm.args,
                                  list(stepmax=1/2, typsize=typsize)))
    popt <- g$estimate
    ## detect error states
    icode <- g$code
    if(!silent && icode > 2) {
      errmess <- nlmcodes[[icode]]
      if(!is.null(errmess)) warning(errmess) else 
      warning("Unrecognised error code ", paste(icode),
              " returned from nlm", call.=FALSE)
    }
    ## return optimised model
    covfunargs[smap] <- popt
    attr(covfunargs, "fitter") <- "ippm"
    fit <- update(fit0, covfunargs=covfunargs, use.internal=TRUE)  
    return(fit)
  }

  ## from help(nlm)
  nlmcodes <- list(c("Relative gradient is close to zero; ",
                     "current iterate is probably solution"),
                   c("Successive iterates are within tolerance; ",
                     "current iterate is probably solution"),
                   c("Last global step failed to locate a point ",
                     "lower than current estimate. ",
                     "Either current estimate is an approximate ",
                     "local minimum of the function ",
                     "or 'steptol' is too small"),
                   "Iteration limit exceeded",
                   c("Maximum step size 'stepmax' ",
                     "exceeded five consecutive times. ",
                     "Either the function is unbounded below, ",
                     "becomes asymptotic to a finite value ",
                     "from above in some direction, ",
                     "or 'stepmax' is too small"))

  ippm
})
