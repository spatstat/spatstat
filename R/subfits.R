#
#
#  $Revision: 1.30 $   $Date: 2014/03/22 01:21:07 $
#
#
subfits.new <- function(object, what="models", verbose=FALSE) {
  stopifnot(inherits(object, "mppm"))

  if(!(what %in% c("models","interactions")))
    stop(paste("Unrecognised option: what=", dQuote(what)))

  if(what == "interactions")
    return(subfits.old(object, what, verbose))
  
  # extract stuff
  announce <- if(verbose) function(x) { cat(x) } else function(x) {} 
    
  announce("Extracting stuff...")
  fitter   <- object$Fit$fitter
  FIT      <- object$Fit$FIT
  coef.FIT <- coef(FIT)
  trend    <- object$trend
  iformula <- object$iformula
  use.gam  <- object$Fit$use.gam
  info     <- object$Info
  npat     <- object$npat
  Inter    <- object$Inter
  interaction <- Inter$interaction
  itags    <- Inter$itags
  Vnamelist <- object$Fit$Vnamelist
  announce("done.\n")

  # determine which interaction(s) are active on each row
  announce("Determining active interactions...")
  active <- active.interactions(object)
  announce("done.\n")

  # exceptions
  if(any(rowSums(active) > 1))
    stop(paste("subfits() is not implemented for models",
               "in which several interpoint interactions",
               "are active on the same point pattern"))
  
  # implied coefficients for each active interaction
  announce("Computing implied coefficients...")
  implcoef <- list()
  for(tag in itags) {
    announce(tag)
    implcoef[[tag]] <- impliedcoefficients(object, tag)
    announce(", ")
  }
  announce("done.\n")

  # Fisher information, if possible
  if(what == "models") {
    announce("Fisher information...")
    fisher   <- vcov(object, what="fisher", err="null")
    varcov   <- try(solve(fisher), silent=TRUE)
    if(inherits(varcov, "try-error"))
      varcov <- NULL
    announce("done.\n")
  }
  
  # Extract data frame 
  announce("Extracting data...")
  datadf   <- object$datadf
  has.design <- info$has.design
  rownames <- object$Info$rownames
  announce("done.\n")

  # set up lists for results 
  models <- rep(list(NULL), npat)
  interactions <- rep(list(NULL), npat)

  # interactions
  announce("Determining interactions...")
  for(i in 1:npat) {
    if(verbose) progressreport(i, npat)
      # Find relevant interaction
    acti <- active[i,]
    nactive <- sum(acti)
    interi <- if(nactive == 0) Poisson else interaction[i, acti, drop=TRUE]
    tagi <- names(interaction)[acti]
      # Find relevant coefficients
    coef.avail  <- coef.FIT
    if(nactive == 1) {
      ic <- implcoef[[tagi]]
      coef.implied <- ic[i, ,drop=TRUE]
      names(coef.implied) <- colnames(ic)
    }
    # overwrite any existing values of coefficients; add new ones.
    coef.avail[names(coef.implied)] <- coef.implied
      # create fitted interaction with these coefficients
    interactions[[i]] <- fii(interi, coef.avail, Vnamelist[[tagi]])
    }
  announce("Done!\n")
  names(interactions) <- rownames

  #
  if(what=="interactions") 
    return(interactions)
  
  # Extract data required to reconstruct complete model fits
  announce("Extracting more data...")
  data  <- object$data
  Y     <- object$Y
  Yname <- info$Yname
  moadf <- object$Fit$moadf
  fmla  <- object$Fit$fmla
  # deal with older formats of mppm
  if(is.null(Yname)) Yname <- info$Xname
  if(is.null(Y)) Y <- data[ , Yname, drop=TRUE]
  # 
  used.cov.names <- info$used.cov.names
  has.covar <- info$has.covar
  if(has.covar) {
    covariates.hf <- data[, used.cov.names, drop=FALSE]
    dfvar <- used.cov.names %in% names(datadf)
  }
  announce("done.\n")

  # Construct template for fake ppm object
  spv <- package_version(versionstring.spatstat())
  fake.version <- list(major=spv$major,
                      minor=spv$minor,
                      release=spv$patchlevel,
                      date="$Date: 2014/03/22 01:21:07 $")
  fake.call <- call("cannot.update", Q=NULL, trend=object$trend,
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
  for(i in 1:npat) {
    if(verbose) progressreport(i, npat)
    Yi <- Y[[i]]
    Wi <- if(is.ppp(Yi)) Yi$window else Yi$data$window
    # assemble relevant covariate images
    covariates <-
      if(has.covar) covariates.hf[i, , drop=TRUE, strip=FALSE] else NULL
    if(has.covar && has.design) 
      ## Convert each data frame covariate value to an image
      covariates[dfvar] <- lapply(covariates[dfvar], as.im, W=Wi)

    # Extract relevant interaction
    finte <- interactions[[i]]
    inte  <- finte$interaction
    if(is.poisson.interact(inte)) inte <- NULL
    Vnames <- finte$Vnames
    if(length(Vnames) == 0) Vnames <- NULL
    
    # Construct fake ppm object
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

    # store in list
    models[[i]] <- fakemodel
  }
  announce("done.\n")
  names(models) <- rownames
  class(models) <- c("listof", class(models))
  return(models)
}

subfits <-
subfits.old <-
  function(object, what="models", verbose=FALSE) {
  stopifnot(inherits(object, "mppm"))

  if(!(what %in% c("models","interactions")))
    stop(paste("Unrecognised option: what=", dQuote(what)))
  
  # extract stuff
  announce <- if(verbose) function(x) { cat(x) } else function(x) {} 
    
  announce("Extracting stuff...")
  FIT      <- object$Fit$FIT
  coef.FIT <- coef(FIT)
  trend    <- object$trend
  iformula <- object$iformula
  use.gam  <- object$Fit$use.gam
  info     <- object$Info
  npat     <- object$npat
  Inter    <- object$Inter
  interaction <- Inter$interaction
  itags    <- Inter$itags
  Vnamelist <- object$Fit$Vnamelist
  announce("done.\n")

  # determine which interaction(s) are active on each row
  announce("Determining active interactions...")
  active <- active.interactions(object)
  announce("done.\n")

  # exceptions
  if(any(rowSums(active) > 1))
    stop(paste("subfits() is not implemented for models",
               "in which several interpoint interactions",
               "are active on the same point pattern"))
  
  # implied coefficients for each active interaction
  announce("Computing implied coefficients...")
  implcoef <- list()
  for(tag in itags) {
    announce(tag)
    implcoef[[tag]] <- impliedcoefficients(object, tag)
    announce(", ")
  }
  announce("done.\n")

  # Fisher information, if possible
  if(what == "models") {
    announce("Fisher information...")
    fisher   <- vcov(object, what="fisher", err="null")
    varcov   <- try(solve(fisher), silent=TRUE)
    if(inherits(varcov, "try-error"))
      varcov <- NULL
    announce("done.\n")
  }
  
  # Extract data frame 
  announce("Extracting data...")
  datadf   <- object$datadf
  has.design <- info$has.design
  rownames <- object$Info$rownames
  announce("done.\n")

  # set up list for results 
  results <- rep(list(NULL), npat)
  
  if(what == "interactions") {
    announce("Determining interactions...")
    for(i in 1:npat) {
      if(verbose) progressreport(i, npat)
      # Find relevant interaction
      acti <- active[i,]
      nactive <- sum(acti)
      interi <- if(nactive == 0) Poisson else interaction[i, acti, drop=TRUE]
      tagi <- names(interaction)[acti]
      # Find relevant coefficients
      coef.avail  <- coef.FIT
      if(nactive == 1) {
        ic <- implcoef[[tagi]]
        coef.implied <- ic[i, ,drop=TRUE]
        names(coef.implied) <- colnames(ic)
      }
      # overwrite any existing values of coefficients; add new ones.
      coef.avail[names(coef.implied)] <- coef.implied
      # create fitted interaction with these coefficients
      results[[i]] <- fii(interi, coef.avail, Vnamelist[[tagi]])
    }
    announce("Done!\n")
    names(results) <- rownames
    return(results)
  }
  
  # Extract data required to reconstruct complete model fits
  announce("Extracting more data...")
  data  <- object$data
  Y     <- object$Y
  Yname <- info$Yname
  # deal with older formats of mppm
  if(is.null(Yname)) Yname <- info$Xname
  if(is.null(Y)) Y <- data[ , Yname, drop=TRUE]
  # 
  used.cov.names <- info$used.cov.names
  has.covar <- info$has.covar
  if(has.covar) {
    covariates.hf <- data[, used.cov.names, drop=FALSE]
    dfvar <- used.cov.names %in% names(datadf)
  }
  announce("done.\n")
  
  ## Loop through point patterns
  announce("Looping through rows...")
  for(i in 1:npat) {
    if(verbose) progressreport(i, npat)
    Yi <- Y[[i]]
    Wi <- if(is.ppp(Yi)) Yi$window else Yi$data$window
    # assemble relevant covariate images
    covariates <-
      if(has.covar) covariates.hf[i, , drop=TRUE, strip=FALSE] else NULL
    if(has.covar && has.design) {
      ## Convert each data frame covariate value to an image
      imrowi <- lapply(covariates[dfvar], as.im, W=Wi)
      
      # Problem: constant covariate leads to singular fit
      # --------------- Hack: ---------------------------
      #  Construct fake data by resampling from possible values
      possible <- function(z) {
        if(is.factor(z))
          factor(levels(z), levels=levels(z))
        else
          unique(z)
      }
      covar.vals <- lapply(as.list(covariates[dfvar, drop=FALSE]), possible)
      scramble <- function(vals, W, Y) {
        W <- as.mask(W)
        npixels <- prod(W$dim)
        nvalues <- length(vals)
        npoints <- Y$n
        # sample the possible values randomly at the non-data pixels
        sampled <- sample(vals, npixels, replace=TRUE)
        Z <- im(sampled, xcol=W$xcol, yrow=W$yrow)
        # repeat the possible values cyclically at the data points
        if(npoints >= 1)
          Z[Y] <- vals[1 + ((1:npoints) %% nvalues)]
        return(Z)
      }
      fake.imrowi <- lapply(covar.vals, scramble, W=Wi, Y=Yi$data)
      # insert fake data into covariates 
      covariates[dfvar] <- fake.imrowi
      # ------------------ end hack ----------------------------
    }

    # Fit ppm to data for case i only
    # using relevant interaction
    acti <- active[i,]
    nactive <- sum(acti)
    if(nactive == 1){
      interi <- interaction[i, acti, drop=TRUE] 
      tagi <- names(interaction)[acti]
      fit <- ppm(Yi, trend, interi, covariates=covariates,
                 use.gam=use.gam,
                 forcefit=TRUE, vnamebase=tagi, vnameprefix=tagi)
    } else {
      fit <- ppm(Yi, trend, Poisson(), covariates=covariates,
                 use.gam=use.gam,
                 forcefit=TRUE)
    }
    
    # now reset the coefficients to those obtained from the full fit
    coefnames.wanted <- names(fit$coef)
    coef.avail  <- coef.FIT
    if(nactive == 1) {
      ic <- implcoef[[tagi]]
      coef.implied <- ic[i, ,drop=TRUE]
      names(coef.implied) <- colnames(ic)
      # overwrite any existing values of coefficients; add new ones.
      coef.avail[names(coef.implied)] <- coef.implied
    }
    if(!all(coefnames.wanted %in% names(coef.avail))) 
      stop("Internal error: some fitted coefficients not accessible")
    fit$theta <- fit$coef <- coef.avail[coefnames.wanted]
    # ... make sure these coefficients will be used in getglmfit, etc ...
    fit$method <- "mppm"
    # ... and replace fake data by true data
    if(has.design) {
      for(nam in names(imrowi)) {
        fit$covariates[[nam]] <- imrowi[[nam]]
        fit$internal$glmdata[[nam]] <- data[i, nam, drop=TRUE]
      }
      # ... and tell glm fit object that it has full rank
      fit$internal$glmfit$rank <- FIT$rank
    }
    # Fisher information and variance-covariance if known
    # Extract submatrices for relevant parameters
    if(!is.null(fisher)) 
      fit$fisher <- fisher[coefnames.wanted, coefnames.wanted, drop=FALSE]
    if(!is.null(varcov))
      fit$varcov <- varcov[coefnames.wanted, coefnames.wanted, drop=FALSE]
    # store in list
    results[[i]] <- fit
  }
  announce("done.\n")
  names(results) <- rownames
  class(results) <- c("listof", class(results))
  return(results)
}

cannot.update <- function(...) {
  stop("This model cannot be updated")
}
