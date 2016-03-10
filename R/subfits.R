#
#
#  $Revision: 1.46 $   $Date: 2016/02/15 12:43:49 $
#
#

subfits.new <- local({

  subfits.new <- function(object, what="models", verbose=FALSE) {
    stopifnot(inherits(object, "mppm"))

    what <- match.arg(what, c("models", "interactions", "basicmodels"))

    if(what == "interactions")
      return(subfits.old(object, what, verbose))
  
    ## extract stuff
    announce <- if(verbose) Announce else Ignore
    
    announce("Extracting stuff...")
    fitter   <- object$Fit$fitter
    FIT      <- object$Fit$FIT
    trend    <- object$trend
#%^!ifdef RANDOMEFFECTS  
    random   <- object$random
#%^!endif  
    info     <- object$Info
    npat     <- object$npat
    Inter    <- object$Inter
    interaction <- Inter$interaction
    itags    <- Inter$itags
    Vnamelist <- object$Fit$Vnamelist
    has.design <- info$has.design
#    has.random <- info$has.random
    announce("done.\n")

    ## fitted parameters
    coefs.full <- coef(object)
    if(is.null(dim(coefs.full))) {
      ## fixed effects model: replicate vector to matrix
      coefs.names <- names(coefs.full)
      coefs.full <- matrix(coefs.full, byrow=TRUE,
                           nrow=npat, ncol=length(coefs.full),
                           dimnames=list(NULL, coefs.names))
    } else {
      ## random/mixed effects model: coerce to matrix
      coefs.names <- colnames(coefs.full)
      coefs.full <- as.matrix(coefs.full)
    }
    
    ## determine which interaction(s) are active on each row
    announce("Determining active interactions...")
    active <- active.interactions(object)
    announce("done.\n")

    ## exceptions
    if(any(rowSums(active) > 1))
      stop(paste("subfits() is not implemented for models",
                 "in which several interpoint interactions",
                 "are active on the same point pattern"))
#%^!ifdef RANDOMEFFECTS  
    if(!is.null(random) && any(variablesinformula(random) %in% itags))
      stop(paste("subfits() is not yet implemented for models",
                 "with random effects that involve",
                 "the interpoint interactions"))
#%^!endif
  
    ## implied coefficients for each active interaction
    announce("Computing implied coefficients...")
    implcoef <- list()
    for(tag in itags) {
      announce(tag)
      implcoef[[tag]] <- impliedcoefficients(object, tag)
      announce(", ")
    }
    announce("done.\n")

    ## Fisher information and vcov
    fisher <- varcov <- NULL
    if(what == "models") {
      announce("Fisher information...")
      fisher   <- vcov(object, what="fisher", err="null")
      varcov   <- try(solve(fisher), silent=TRUE)
      if(inherits(varcov, "try-error"))
        varcov <- NULL
      announce("done.\n")
    } 
  
    ## Extract data frame 
    announce("Extracting data...")
    datadf   <- object$datadf
    rownames <- object$Info$rownames
    announce("done.\n")

    ## set up lists for results 
    models <- rep(list(NULL), npat)
    interactions <- rep(list(NULL), npat)
    
    ## interactions
    announce("Determining interactions...")
    pstate <- list()
    for(i in 1:npat) {
      if(verbose) pstate <- progressreport(i, npat, state=pstate)
      ## Find relevant interaction
      acti <- active[i,]
      nactive <- sum(acti)
      interi <- if(nactive == 0) Poisson() else interaction[i, acti, drop=TRUE]
      tagi <- names(interaction)[acti]
      ## Find relevant coefficients
      coefs.avail  <- coefs.full[i,]
      names(coefs.avail) <- coefs.names
      if(nactive == 1) {
        ic <- implcoef[[tagi]]
        coefs.implied <- ic[i, ,drop=TRUE]
        names(coefs.implied) <- colnames(ic)
        ## overwrite any existing values of coefficients; add new ones.
        coefs.avail[names(coefs.implied)] <- coefs.implied
      }
      ## create fitted interaction with these coefficients
      vni <- if(nactive > 0) Vnamelist[[tagi]] else character(0)
      interactions[[i]] <- fii(interi, coefs.avail, vni)
    }
    announce("Done!\n")
    names(interactions) <- rownames

    ##
    if(what=="interactions") 
      return(interactions)
  
    ## Extract data required to reconstruct complete model fits
    announce("Extracting more data...")
    data  <- object$data
    Y     <- object$Y
    Yname <- info$Yname
    moadf <- object$Fit$moadf
    fmla  <- object$Fit$fmla
    ## deal with older formats of mppm
    if(is.null(Yname)) Yname <- info$Xname
    if(is.null(Y)) Y <- data[ , Yname, drop=TRUE]
    ## 
    used.cov.names <- info$used.cov.names
    has.covar <- info$has.covar
    if(has.covar) {
      covariates.hf <- data[, used.cov.names, drop=FALSE]
      dfvar <- used.cov.names %in% names(datadf)
    }
    announce("done.\n")

    ## Construct template for fake ppm object
    spv <- package_version(versionstring.spatstat())
    fake.version <- list(major=spv$major,
                         minor=spv$minor,
                         release=spv$patchlevel,
                         date="$Date: 2016/02/15 12:43:49 $")
    fake.call <- call("cannot.update", Q=NULL, trend=trend,
                      interaction=NULL, covariates=NULL,
                      correction=object$Info$correction,
                      rbord     = object$Info$rbord)
    fakemodel <- list(
                    method       = "mpl",
                    fitter       = fitter,
                    coef         = coef(object),
                    trend        = object$trend,
                    interaction  = NULL,
                    fitin        = NULL,
                    Q            = NULL,
                    maxlogpl     = NA,
                    internal     = list(glmfit = FIT,
                                        glmdata  = NULL,
                                        Vnames   = NULL,
                                        fmla     = fmla,
                                        computed = list()),
                    covariates   = NULL,
                    correction   = object$Info$correction,
                    rbord        = object$Info$rbord,
                    version      = fake.version,
                    problems     = list(),
                    fisher       = fisher,
                    varcov       = varcov,
                    call         = fake.call,
                    callstring   = "cannot.update()",
                    fake         = TRUE)
    class(fakemodel) <- "ppm"

    ## Loop through point patterns
    announce("Generating models for each row...")
    pstate <- list()
    for(i in 1:npat) {
      if(verbose) pstate <- progressreport(i, npat, state=pstate)
      Yi <- Y[[i]]
      Wi <- if(is.ppp(Yi)) Yi$window else Yi$data$window
      ## assemble relevant covariate images
      covariates <-
        if(has.covar) covariates.hf[i, , drop=TRUE, strip=FALSE] else NULL
      if(has.covar && has.design) 
        ## Convert each data frame covariate value to an image
        covariates[dfvar] <- lapply(covariates[dfvar], as.im, W=Wi)

      ## Extract relevant interaction
      finte <- interactions[[i]]
      inte  <- finte$interaction
      if(is.poisson.interact(inte)) inte <- NULL
      Vnames <- finte$Vnames
      if(length(Vnames) == 0) Vnames <- NULL
    
      ## Construct fake ppm object
      fakemodel$interaction <- inte
      fakemodel$fitin       <- finte
      fakemodel$Q           <- Yi
      fakemodel$covariates  <- covariates
      fakemodel$internal$glmdata <- moadf[moadf$id == i, ]
      fakemodel$internal$Vnames  <- Vnames

      fake.call$Q <- Yi
      fake.call$covariates <- covariates
      fakemodel$call <- fake.call
      fakemodel$callstring <- short.deparse(fake.call)
      
      ## store in list
      models[[i]] <- fakemodel
    }
    announce("done.\n")
    names(models) <- rownames
    models <- as.anylist(models)
    return(models)
  }

  Announce <- function(...) cat(...)

  Ignore <- function(...) { NULL }

  subfits.new
})



## /////////////////////////////////////////////////////

subfits <-
subfits.old <- local({
    
  subfits.old <- function(object, what="models", verbose=FALSE) {
    stopifnot(inherits(object, "mppm"))
    what <- match.arg(what, c("models","interactions", "basicmodels"))
    ## extract stuff
    announce <- if(verbose) Announce else Ignore
    
    announce("Extracting stuff...")
#    FIT      <- object$Fit$FIT
    trend    <- object$trend
    random   <- object$random
    use.gam  <- object$Fit$use.gam
    info     <- object$Info
    npat     <- object$npat
    Inter    <- object$Inter
    interaction <- Inter$interaction
    itags    <- Inter$itags
    Vnamelist <- object$Fit$Vnamelist
    has.design <- info$has.design
    has.random <- info$has.random
    announce("done.\n")

    ## fitted parameters
    coefs.full <- coef(object)
    if(is.null(dim(coefs.full))) {
      ## fixed effects model: replicate vector to matrix
      coefs.names <- names(coefs.full)
      coefs.full <- matrix(coefs.full, byrow=TRUE,
                           nrow=npat, ncol=length(coefs.full),
                           dimnames=list(NULL, coefs.names))
    } else {
      ## random/mixed effects model: coerce to matrix
      coefs.names <- colnames(coefs.full)
      coefs.full <- as.matrix(coefs.full)
    }
  
    ## determine which interaction(s) are active on each row
    announce("Determining active interactions...")
    active <- active.interactions(object)
    announce("done.\n")

    ## exceptions
    if(any(rowSums(active) > 1))
      stop(paste("subfits() is not implemented for models",
                 "in which several interpoint interactions",
                 "are active on the same point pattern"))
#%^!ifdef RANDOMEFFECTS  
    if(!is.null(random) && any(variablesinformula(random) %in% itags))
      stop(paste("subfits() is not yet implemented for models",
                 "with random effects that involve",
                 "the interpoint interactions"))
#%^!endif
  
    ## implied coefficients for each active interaction
    announce("Computing implied coefficients...")
    implcoef <- list()
    for(tag in itags) {
      announce(tag)
      implcoef[[tag]] <- impliedcoefficients(object, tag)
      announce(", ")
    }
    announce("done.\n")

    ## Fisher information and vcov
    fisher <- varcov <- NULL
    if(what == "models") {
      announce("Fisher information...")
      fisher   <- vcov(object, what="fisher", err="null")
      varcov   <- try(solve(fisher), silent=TRUE)
      if(inherits(varcov, "try-error"))
        varcov <- NULL
      announce("done.\n")
    }
  
    ## Extract data frame 
    announce("Extracting data...")
    datadf   <- object$datadf
    rownames <- object$Info$rownames
    announce("done.\n")

    ## set up list for results 
    results <- rep(list(NULL), npat)
  
    if(what == "interactions") {
      announce("Determining interactions...")
      pstate <- list()
      for(i in 1:npat) {
        if(verbose) pstate <- progressreport(i, npat, state=pstate)
        ## Find relevant interaction
        acti <- active[i,]
        nactive <- sum(acti)
        interi <- if(nactive == 0) Poisson() else
                  interaction[i, acti, drop=TRUE]
        tagi <- names(interaction)[acti]
        ## Find relevant coefficients
        coefs.avail  <- coefs.full[i,]
        names(coefs.avail) <- coefs.names
        if(nactive == 1) {
          ic <- implcoef[[tagi]]
          coefs.implied <- ic[i, ,drop=TRUE]
          names(coefs.implied) <- colnames(ic)
          ## overwrite any existing values of coefficients; add new ones.
          coefs.avail[names(coefs.implied)] <- coefs.implied
        }
        ## create fitted interaction with these coefficients
        vni <- if(nactive > 0) Vnamelist[[tagi]] else character(0)
        results[[i]] <- fii(interi, coefs.avail, vni)
      }
      announce("Done!\n")
      names(results) <- rownames
      return(results)
    }
  
    ## Extract data required to reconstruct complete model fits
    announce("Extracting more data...")
    data  <- object$data
    Y     <- object$Y
    Yname <- info$Yname
    ## deal with older formats of mppm
    if(is.null(Yname)) Yname <- info$Xname
    if(is.null(Y)) Y <- data[ , Yname, drop=TRUE]
    ##
    used.cov.names <- info$used.cov.names
    has.covar <- info$has.covar
    if(has.covar) {
      covariates.hf <- data[, used.cov.names, drop=FALSE]
      dfvar <- used.cov.names %in% names(datadf)
    }
    announce("done.\n")
  
    ## Loop through point patterns
    announce("Looping through rows...")
    pstate <- list()
    for(i in 1:npat) {
      if(verbose) pstate <- progressreport(i, npat, state=pstate)
      Yi <- Y[[i]]
      Wi <- if(is.ppp(Yi)) Yi$window else Yi$data$window
      ## assemble relevant covariate images
      covariates <-
        if(has.covar) covariates.hf[i, , drop=TRUE, strip=FALSE] else NULL
      if(has.covar && has.design) {
        ## Convert each data frame covariate value to an image
        imrowi <- lapply(covariates[dfvar], as.im, W=Wi)
        ## Problem: constant covariate leads to singular fit
        ## --------------- Hack: ---------------------------
        ##  Construct fake data by resampling from possible values
        covar.vals <- lapply(as.list(covariates[dfvar, drop=FALSE]), possible)
        fake.imrowi <- lapply(covar.vals, scramble, W=Wi, Y=Yi$data)
        ## insert fake data into covariates 
        covariates[dfvar] <- fake.imrowi
        ## ------------------ end hack ----------------------------
      }
      ## Fit ppm to data for case i only
      ## using relevant interaction
      acti <- active[i,]
      nactive <- sum(acti)
      if(nactive == 1){
        interi <- interaction[i, acti, drop=TRUE] 
        tagi <- names(interaction)[acti]
        fiti <- PiPiM(Yi, trend, interi, covariates=covariates,
                      allcovar=has.random,
                      use.gam=use.gam,
                      vnamebase=tagi, vnameprefix=tagi)
      } else {
        fiti <- PiPiM(Yi, trend, Poisson(), covariates=covariates,
                      allcovar=has.random,
                      use.gam=use.gam)
      }
      ## fiti determines which coefficients are required
      coefi.fitted <- fiti$coef
      coefnames.wanted <- names(coefi.fitted)
      ## take the required coefficients from the full mppm fit
      coefs.avail  <- coefs.full[i,]
      names(coefs.avail) <- coefs.names
      if(nactive == 1) {
        ic <- implcoef[[tagi]]
        coefs.implied <- ic[i, ,drop=TRUE]
        names(coefs.implied) <- colnames(ic)
        ## overwrite any existing values of coefficients; add new ones.
        coefs.avail[names(coefs.implied)] <- coefs.implied
      }
      if(!all(coefnames.wanted %in% names(coefs.avail))) 
        stop("Internal error: some fitted coefficients not accessible")
      coefi.new <- coefs.avail[coefnames.wanted]
      ## reset coefficients
      fiti$coef.orig <- coefi.fitted ## (detected by summary.ppm, predict.ppm)
      fiti$theta <- fiti$coef <- coefi.new
      fiti$method <- "mppm"
      ## ... and replace fake data by true data
      if(has.design) {
        for(nam in names(imrowi)) {
          fiti$covariates[[nam]] <- imrowi[[nam]]
          fiti$internal$glmdata[[nam]] <- data[i, nam, drop=TRUE]
        }
      }
      ## Adjust rank of glm fit object
#      fiti$internal$glmfit$rank <- FIT$rank 
      fiti$internal$glmfit$rank <- sum(is.finite(fiti$coef))
      ## Fisher information and variance-covariance if known
      ## Extract submatrices for relevant parameters
      if(!is.null(fisher)) 
        fiti$fisher <- fisher[coefnames.wanted, coefnames.wanted, drop=FALSE]
      if(!is.null(varcov))
        fiti$varcov <- varcov[coefnames.wanted, coefnames.wanted, drop=FALSE]
      ## store in list
      results[[i]] <- fiti
    }
    announce("done.\n")
    names(results) <- rownames
    results <- as.anylist(results)
    return(results)
  }

  PiPiM <- function(Y, trend, inter, covariates, ...,
                    allcovar=FALSE, use.gam=FALSE,
                    vnamebase=c("Interaction", "Interact."),
                    vnameprefix=NULL) {
    # This ensures that the model is fitted in a unique environment
    # so that it can be updated later.
    force(Y)
    force(trend)
    force(inter)
    force(covariates)
    force(allcovar)
    force(use.gam)
    force(vnamebase)
    force(vnameprefix)
    feet <- ppm(Y, trend, inter, covariates=covariates,
                allcovar=allcovar, use.gam=use.gam,
                forcefit=TRUE, vnamebase=vnamebase, vnameprefix=vnameprefix)
    return(feet)
  }
  
  possible <- function(z) {
    if(!is.factor(z)) unique(z) else factor(levels(z), levels=levels(z))
  }
  
  scramble <- function(vals, W, Y) {
    W <- as.mask(W)
    npixels <- prod(W$dim)
    nvalues <- length(vals)
    npts <- npoints(Y)
    ## sample the possible values randomly at the non-data pixels
    sampled <- sample(vals, npixels, replace=TRUE)
    Z <- im(sampled, xcol=W$xcol, yrow=W$yrow)
    ## repeat the possible values cyclically at the data points
    if(npts >= 1)
      Z[Y] <- vals[1 + ((1:npts) %% nvalues)]
    return(Z)
  }

  Announce <- function(...) cat(...)

  Ignore <- function(...) { NULL }

  subfits.old
})

cannot.update <- function(...) {
  stop("This model cannot be updated")
}
