#
# ippm.R
#
#   $Revision: 2.7 $   $Date: 2011/12/03 11:37:03 $
#
# Fisher scoring algorithm for irregular parameters in ppm trend
#

ippm <- function(...,
                   iScore=NULL, 
                   start=list(),
                   covfunargs=start,
                   maxiter=20, tol=1e-4, progress=TRUE, stepfactor=1,
                   dbug=FALSE) {
  # validate
  if(!is.list(start) || length(start) == 0)
    stop("start should be a list of initial values for irregular parameters")
  if(!is.list(iScore) || length(iScore) != length(start))
    stop("iScore should be a list of the same length as start")
  stopifnot(identical(names(iScore), names(start)))
  if(!all(unlist(lapply(iScore, is.function))))
    stop("iScore should be a list of functions")
  #
  smap <- match(names(start), names(covfunargs))
  if(any(is.na(smap)))
      stop("variables in start should be a subset of variables in covfunargs")
  covfunargs[smap] <- start
  #
  # fit the initial model and extract information
  fit0 <- ppm(..., covfunargs=covfunargs)
  lpl0 <- fit0$maxlogpl
  p <- length(coef(fit0))
  # examine covariates and trend
  covariates <- fit0$covariates
  isfun <- unlist(lapply(covariates, is.function))
  covfuns <- covariates[isfun]
  # determine which covariates depend on which irregular parameters
  pnames <- names(start)
  hasarg <- function(f,a) { a %in% names(formals(f)) }
  depmat <- matrix(FALSE, nrow=length(covfuns), ncol=length(pnames))
  rownames(depmat) <- names(covfuns)
  colnames(depmat) <- pnames
  for(j in 1:length(pnames))
    depmat[,j] <- unlist(lapply(covfuns, hasarg, pnames[j]))
  # find covariates that depend on ANY irregular parameter 
  depvar <- rownames(depmat)[apply(depmat, 1, any)]
  # check that these covariates appear only in offset terms
  covnames.fitted <- model.covariates(fit0, fitted=TRUE,  offset=FALSE)
  if(any(uhoh <- depvar %in% covnames.fitted))
    stop(paste(ngettext(sum(uhoh), "The covariate", "The covariates"),
               commasep(sQuote(depvar[uhoh])),
               "should appear only in offset terms"))
  # check that every irregular parameter to be updated appears somewhere 
  covnames.offset <- model.covariates(fit0, fitted=FALSE,  offset=TRUE)
  usearg <- apply(depmat[covnames.offset, , drop=FALSE], 2, any)
  if(!all(usearg)) {
    nbad <- sum(!usearg)
    warning(paste("Cannot maximise over the irregular",
               ngettext(nbad, "parameter", "parameters"),
               commasep(sQuote(names(usearg)[!usearg])),
               ngettext(nbad, "because it is", "because they are"),
               "not used in any term of the model"))
    # restrict 
    start <- start[usearg]
    iScore <- iScore[usearg]
    pnames <- names(start)
  }
  # ready
  iter <- 0
  param <- start
  pvec <- as.numeric(param)
  ndigits <- max(0, -ceiling(log10(tol))) + 1
  # go
  for(phase in 1:2) {
    if(progress) cat(paste("Phase", phase, "\n"))
    maxit <- if(phase == 1) 2 else maxiter
    for(iter in 0:maxit) {
    # fit model with current irregular parameters
      covfunargs[smap] <- param
      fit <- ppm(..., covfunargs=covfunargs)
      lpl <- logLik(fit, warn=FALSE)
      if(progress) {
        co <- coef(fit)
        cat(paste(paren(iter, "["),
                  paste(paste(pnames, "=", round(pvec, digits=ndigits)),
                        collapse=", "),
                  "->",
                  paste(paste(names(co), "=", round(co, digits=ndigits)),
                        collapse=", "),
                  "; logPL=", signif(lpl, 8), "\n"))
      }
      # compute model matrix and inverse fisher information
      stuff <- ppm.influence(fit, what="derivatives",
                             iScore=iScore,
                             iArgs=param)
      score <- stuff$deriv$score
      vc    <- stuff$deriv$vc
      fush  <- stuff$deriv$fush
      if(dbug) {
        cat("\nscore=\n")
        print(score)
        cat("\nvc=\n")
        print(vc)
        cat("\nfush=\n")
        print(fush)
      }
      if(phase == 1) {
        # Fisher scoring on partial matrix 
        subscore <- score[, -(1:p), drop=FALSE]
        subfush <- fush[-(1:p), -(1:p)]
        stepvec <- as.numeric(solve(subfush) %*% t(subscore))
      } else {
        # Fisher scoring step on full parameter vector 
        stepvec <- as.numeric(vc %*% t(score))[ -(1:p)]
      }
      if(dbug) {
        cat("\nstep=\n")
        print(stepvec)
      }
      # update parameters or exit
      if(iter > 0 && all(abs(stepvec) < tol))
        break
      pvec <- pvec + stepfactor * stepvec
      param <- as.list(pvec)
      names(param) <- pnames
    }
  }
  if(iter == maxiter)
    warning(paste("Maximum number of iterations", paren(maxiter),
                  "reached without convergence\n"))
  return(fit)
}
