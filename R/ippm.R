#
# ippm.R
#
#   $Revision: 2.21 $   $Date: 2017/02/07 07:47:20 $
#
# Fisher scoring algorithm for irregular parameters in ppm trend
#

ippm <- local({

  chucknames <- c("iScore", "start", "nlm.args", "silent", "warn.unused")
  
  ippm <- function(Q, ...,
                   iScore=NULL, 
                   start=list(),
                   covfunargs=start,
                   nlm.args=list(),
                   silent=FALSE,
                   warn.unused=TRUE) {
    ## remember call
    cl <- match.call()
    callframe <- parent.frame()
    callstring <- short.deparse(sys.call())
    ##
    ppmcall <- cl[!(names(cl) %in% chucknames)]
    ppmcall[[1L]] <- as.name('ppm')
    ## validate
    if(!is.list(start))
      stop("start should be a list of initial values for irregular parameters")
    if(length(start) == 0) {
      ppmcall <- ppmcall[names(ppmcall) != "covfunargs"]
      return(eval(ppmcall, callframe))
    }
    if(!is.null(iScore)) {
      if(!is.list(iScore) || length(iScore) != length(start))
        stop("iScore should be a list of the same length as start")
      stopifnot(identical(names(iScore), names(start)))
      if(!all(unlist(lapply(iScore, is.function))))
        stop("iScore should be a list of functions")
    }
    ##
    smap <- match(names(start), names(covfunargs))
    if(anyNA(smap))
      stop("variables in start should be a subset of variables in covfunargs")
    covfunargs[smap] <- start
    ## fit the initial model and extract information
    ppmcall$covfunargs <- covfunargs
    fit0 <- eval(ppmcall, callframe)
#    lpl0 <- fit0$maxlogpl
#    p <- length(coef(fit0))
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
    depvar <- rownames(depmat)[apply(depmat, 1L, any)]
    ## check that these covariates appear only in offset terms
    covnames.fitted <- model.covariates(fit0, fitted=TRUE,  offset=FALSE)
    if(any(uhoh <- depvar %in% covnames.fitted))
      stop(paste(ngettext(sum(uhoh), "The covariate", "The covariates"),
                 commasep(sQuote(depvar[uhoh])),
                 "should appear only in offset terms"))
    ## check that every irregular parameter to be updated appears somewhere 
    cov.names.offset <- model.covariates(fit0, fitted=FALSE,  offset=TRUE)
    covfun.names.offset <- intersect(cov.names.offset, names(covfuns))
    usearg <- apply(depmat[covfun.names.offset, , drop=FALSE], 2L, any)
    if(!all(usearg)) {
      if(warn.unused) {
        nbad <- sum(!usearg)
        warning(paste("Cannot maximise over the irregular",
                      ngettext(nbad, "parameter", "parameters"),
                      commasep(sQuote(names(usearg)[!usearg])),
                      ngettext(nbad, "because it is", "because they are"),
                      "not used in any term of the model"))
      }
      ## restrict 
      start <- start[usearg]
      if(!is.null(iScore)) iScore <- iScore[usearg]
      pnames <- names(start)
    }
    if(length(start) == 0) {
      ppmcall <- ppmcall[names(ppmcall) != "covfunargs"]
      return(eval(ppmcall, callframe))
    }
    ## parameters for objective function
    fdata <- list(fit0=fit0,
                  nreg=length(coef(fit0)),
                  covfunargs=covfunargs,
                  smap=smap,
                  pnames=pnames,
                  iScore=iScore)
    ## minimise objective
    startvec <- unlist(start)
    typsize <- abs(startvec)
    typsize <- pmax(typsize, min(typsize[typsize > 0]))
    g <- do.call(nlm,
                 resolve.defaults(list(f=objectivefun,
                                       p=startvec,
                                       thedata=fdata),
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
    attr(covfunargs, "free") <- names(start)
    fit <- update(fit0, covfunargs=covfunargs, use.internal=TRUE)
    fit$dispatched <- fit[c("call", "callstring", "callframe")]
    fit$call <- cl
    fit$callstring <- callstring
    fit$callframe <- callframe
    class(fit) <- c("ippm", class(fit))
    return(fit)
  }

  ## define objective function
  objectivefun <- function(param, thedata) {
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

  ## from help(nlm)
  nlmcodes <- list(c("Relative gradient is close to zero; ",
                     "current iterate is probably solution"),
                   c("Successive iterates are within tolerance; ",
                     "current iterate is probably solution"),
                   c("Last global step failed to locate a point ",
                     "lower than current estimate. ",
                     "Either current estimate is an approximate ",
                     "local minimum of the function ",
                     "or 'steptol' is too large"),
                   "Iteration limit exceeded",
                   c("Maximum step size 'stepmax' ",
                     "exceeded five consecutive times. ",
                     "Either the function is unbounded below, ",
                     "becomes asymptotic to a finite value ",
                     "from above in some direction, ",
                     "or 'stepmax' is too small"))

  ippm
})


update.ippm <- local({

  newformula <- function(old, change, eold, enew) {
    old <- eval(old, eold)
    change <- eval(change, enew)
    old <- as.formula(old, env=eold)
    change <- as.formula(change, env=enew)
    update.formula(old, change)
  }

  update.ippm <- function(object, ..., envir=environment(terms(object))) {
#    call <- match.call()
    new.call <- old.call <- object$call
    old.callframe <- object$callframe
    Qold <- eval(old.call$Q, as.list(envir), enclos=old.callframe)
    argh <- list(...)
    if(any(isfmla <- unlist(lapply(argh, inherits, what="formula")))) {
      if(sum(isfmla) > 1)
        stop("Syntax not understood: several arguments are formulas")
      i <- min(which(isfmla))
      new.fmla <- argh[[i]]
      argh <- argh[-i]
      if(inherits(Qold, "formula")) {
        ## formula will replace 'Q'
        if(is.null(lhs.of.formula(new.fmla))) {
          f <- (. ~ x)
          f[[3L]] <- new.fmla[[2L]]
          new.fmla <- f
        }
        new.call$Q <- newformula(Qold, new.fmla, old.callframe, envir)
      } else if(inherits(Qold, c("ppp", "quad"))) {
        ## formula will replace 'trend' and may replace 'Q'
        new.fmla <- newformula(formula(object), new.fmla, old.callframe, envir)
        if(!is.null(lhs <- lhs.of.formula(new.fmla))) {
          newQ <- eval(eval(substitute(substitute(l, list("."=Q)),
                                       list(l=lhs,
                                            Q=Qold))),
                       envir=as.list(envir), enclos=old.callframe)
          new.call$Q <- newQ
        }
        new.fmla <- rhs.of.formula(new.fmla)
        if("trend" %in% names(old.call)) {
          new.call$trend <- new.fmla
        } else {
          ## find which argument in the original call was a formula
          wasfmla <- sapply(old.call, formulaic,
                            envir=as.list(envir),
                            enclos=old.callframe)
          if(any(wasfmla)) {
            new.call[[min(which(wasfmla))]] <- new.fmla
          } else {
            new.call$trend <- new.fmla
          }
        }
      }
    }
    ## silence the warnings about unused covfunargs (unless overruled)
    new.call$warn.unused <- FALSE
    ## other arguments
    if(length(argh) > 0) {
      nama <- names(argh)
      named <- if(is.null(nama)) rep(FALSE, length(argh)) else nzchar(nama)
      if(any(named))
        new.call[nama[named]] <- argh[named]
      if(any(!named))
        new.call[length(new.call) + 1:sum(!named)] <- argh[!named]
    }
    result <- eval(new.call, as.list(envir), enclos=old.callframe)
    return(result)
  }

  formulaic <- function(z, envir, enclos) {
    u <- try(eval(z, envir, enclos))
    return(inherits(u, "formula"))
  }

  update.ippm
})

