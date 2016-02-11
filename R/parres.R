#
# parres.R
#
# code to plot transformation diagnostic
#
#   $Revision: 1.5 $  $Date: 2016/02/11 10:17:12 $
#

parres <- function(model, covariate, ...,
                   smooth.effect=FALSE, subregion=NULL,
                   bw="nrd0", adjust=1, from=NULL,to=NULL, n=512,
                   bw.input = c("points", "quad"),
                   bw.restrict = FALSE,
                   covname) {  

  modelname <- deparse(substitute(model))
  if(missing(covname)) 
    covname <- sensiblevarname(deparse(substitute(covariate)), "X")
  callstring <- paste(deparse(sys.call()), collapse = "")

  if(!is.null(subregion)) 
    stopifnot(is.owin(subregion))
  
  if(is.null(adjust)) adjust <- 1

  bw.input <- match.arg(bw.input)
  
  # validate model
  stopifnot(is.ppm(model))
  modelcall <- model$callstring
  if(is.null(modelcall))
    modelcall <- model$call
  if(is.null(getglmfit(model)))
    model <- update(model, forcefit=TRUE)
  
  # extract spatial locations
  Q <- quad.ppm(model)
#  datapoints <- Q$data
  quadpoints <- union.quad(Q)
  Z <- is.data(Q)
  wts <- w.quad(Q)
  nQ <- npoints(quadpoints)
  # fitted intensity
  lam <- fitted(model, type="trend")
  # subset of quadrature points used to fit model
  subQset <- getglmsubset(model)
  if(is.null(subQset)) subQset <- rep.int(TRUE, nQ)
  # restriction to subregion
  insubregion <- if(!is.null(subregion)) {
    inside.owin(quadpoints, w=subregion)
  } else rep.int(TRUE, nQ)

  ################################################################
  # Inverse lambda residuals

  rx <- residuals(model, type="inverse")
  resid <- with(rx, "increment")

  #################################################################
  # identify the covariate
  #
  if(length(covariate) == 0)
    stop("No covariate specified")

  covtype <- "unknown"

  if(!is.character(covariate)) {
    # Covariate is some kind of data, treated as external covariate
    covtype <- "external"
    beta <- 0
    covvalues <- evalCovariate(covariate, quadpoints)
  } else {
    # Argument is name of covariate
    covname <- covariate
    if(length(covname) > 1)
      stop("Must specify only one covariate")
    # 'original covariates'
    orig.covars <- variablesinformula(formula(model))
    # 'canonical covariates'
    canon.covars <- names(coef(model))
    # offsets
    offset.covars <- offsetsinformula(formula(model))
    # 
    if(covname %in% orig.covars) {
      # one of the original covariates
      covtype <- "original"
      covvalues <- evalCovariate(findCovariate(covname, model), quadpoints)
    } else if(covname %in% canon.covars) {
      # one of the canonical covariates
      covtype <- "canonical"
      mm <- model.matrix(model)
      covvalues <- mm[, covname]
      ## extract the corresponding coefficient
      beta <- coef(model)[[covname]]
    } else if(covname %in% offset.covars) {
      # an offset term only
      covtype <- "offset"
      mf <- model.frame(model, subset=rep.int(TRUE, n.quad(Q)))
      if(!(covname %in% colnames(mf)))
        stop(paste("Internal error: offset term", covname,
                   "not found in model frame"))
      covvalues <- mf[, covname]
      ## fixed coefficient (not an estimated parameter)
      beta <- 1
    } else{
      # must be an external covariate (i.e. not used in fitted model)
      covtype <- "external"
      beta <- 0
      covvalues <- evalCovariate(findCovariate(covname, model), quadpoints)
    }
  }
  # validate covvalues
  #
  if(is.null(covvalues))
    stop("Unable to extract covariate values")
  if(length(covvalues) != npoints(quadpoints))
    stop(paste("Internal error: number of covariate values =",
               length(covvalues), "!=", npoints(quadpoints),
               "= number of quadrature points"))
  vtype <- typeof(covvalues)
  switch(vtype,
         real=,
         double = { },
         integer = {
           warning("Covariate is integer-valued")
         },
         stop(paste("Cannot handle covariate of type", sQuote(vtype))))
  
  #################################################################
  # Compute covariate effect

  if(covtype != "original") {
    effect <- beta * covvalues
    mediator <- covtype
    effectfundata <- list(beta=beta)
    effectFun <- function(x) { (effectfundata$beta) * x }
    isoffset <- (covtype == "offset")
    names(isoffset) <- covname
  } else {
    ## `original' covariate (passed as argument to ppm)
    ## may determine one or more canonical covariates and/or offsets
    origcovdf <- getppmOriginalCovariates(model)[insubregion, , drop=FALSE]
    isconstant <- lapply(origcovdf,
                         function(z) { length(unique(z)) == 1 })
    ##
    ## Initialise
    termnames <- character(0)
    termbetas <- numeric(0)
    isoffset <- logical(0)
    mediator <- character(0)
    effect <- 0
    effectFun <- function(x) { effectFun.can(x) + effectFun.off(x) }
    effectFun.can <- effectFun.off <- function(x) { 0 * x }
    ## Identify relevant canonical covariates
    dmat <- model.depends(model)
    if(!(covname %in% colnames(dmat)))
      stop("Internal error: cannot match covariate names")
    othercov <- (colnames(dmat) != covname)
    relevant <- dmat[, covname]
    if(any(relevant)) {
      # original covariate determines one or more canonical covariates
      mediator <- "canonical"
      # check whether covariate is separable
      if(any(conflict <- dmat[relevant, othercov, drop=FALSE])) {
        ## identify entangled covariates
        entangled <- colnames(conflict)[apply(conflict, 2, any)]
        ## not problematic if constant
        ok <- unlist(isconstant[entangled])
        conflict[ , ok] <- FALSE
        ## re-test
        if(any(conflict)) {
          conflictterms <- apply(conflict, 1, any)
          conflictcovs  <- apply(conflict, 2, any)
          stop(paste("The covariate", sQuote(covname),
                     "cannot be separated from the",
                     ngettext(sum(conflictcovs), "covariate", "covariates"),
                     commasep(sQuote(colnames(conflict)[conflictcovs])),
                     "in the model",
                     ngettext(sum(conflictterms), "term", "terms"),
                     commasep(sQuote(rownames(conflict)[conflictterms]))
                     ))
        }
      }
      # 
      termnames <- rownames(dmat)[relevant]
      isoffset <- rep.int(FALSE, length(termnames))
      names(isoffset) <- termnames
      # Extract relevant canonical covariates
      mm <-  model.matrix(model)
      termvalues <- mm[, relevant, drop=FALSE]
      # extract corresponding coefficients
      termbetas <- coef(model)[relevant]
      # evaluate model effect
      effect <- as.numeric(termvalues %*% termbetas)
      # check length
      if(length(effect) != npoints(quadpoints))
        stop(paste("Internal error: number of values of fitted effect =",
                   length(effect), "!=", npoints(quadpoints),
                   "= number of quadrature points"))
      # Trap loglinear case
      if(length(termnames) == 1 && identical(termnames, covname)) {
        covtype <- "canonical"
        beta <- termbetas
      }
      # construct the corresponding function
      gd <- getglmdata(model)
      goodrow <- min(which(complete.cases(gd)))
      defaultdata <- gd[goodrow, , drop=FALSE]
      effectfundata.can <- list(covname=covname,
                            fmla = formula(model),
                            termbetas = termbetas,
                            defaultdata = defaultdata,
                            relevant = relevant,
                            termnames = termnames)
      effectFun.can <- function(x) {
        d <- effectfundata.can
        # replicate default data to correct length
        df <- as.data.frame(lapply(d$defaultdata, rep, length(x)))
        # overwrite value of covariate with new data
        df[,covname] <- x
        # construct model matrix 
        m <- model.matrix(d$fmla, df)
        # check it conforms to expected structure
        if(!identical(colnames(m)[d$relevant], d$termnames))
          stop("Internal error: mismatch in term names in effectFun")
        me <- m[, d$relevant, drop=FALSE]
        y <- me %*% as.matrix(d$termbetas, ncol=1) 
        return(y)
      }
    }
    if(!is.null(offmat <- attr(dmat, "offset")) &&
       any(relevant <- offmat[, covname])) {
      # covariate appears in a model offset term
      mediator <- c(mediator, "offset")
      # check whether covariate is separable
      if(any(conflict<- offmat[relevant, othercov, drop=FALSE])) {
        ## identify entangled covariates
        entangled <- colnames(conflict)[apply(conflict, 2, any)]
        ## not problematic if constant
        ok <- unlist(isconstant[entangled])
        conflict[ , ok] <- FALSE
        ## re-test
        if(any(conflict)) {
          conflictterms <- apply(conflict, 1, any)
          conflictcovs  <- apply(conflict, 2, any)
          stop(paste("The covariate", sQuote(covname),
                     "cannot be separated from the",
                     ngettext(sum(conflictcovs), "covariate", "covariates"),
                     commasep(sQuote(colnames(conflict)[conflictcovs])),
                     "in the model",
                     ngettext(sum(conflictterms), "term", "terms"),
                     commasep(sQuote(rownames(conflict)[conflictterms]))
                     ))
        }
      }
      # collect information about relevant offset 
      offnames <- rownames(offmat)[relevant]
      termnames <- c(termnames, offnames)
      noff <- length(offnames)
      termbetas <- c(termbetas, rep.int(1, noff))
      isoffset  <- c(isoffset, rep.int(TRUE, noff))
      names(termbetas) <- names(isoffset) <- termnames
      # extract values of relevant offset 
      mf <- model.frame(model, subset=rep.int(TRUE, n.quad(Q)))
      if(any(nbg <- !(offnames %in% colnames(mf))))
        stop(paste("Internal error:",
                   ngettext(sum(nbg), "offset term", "offset terms"),
                   offnames[nbg],
                   "not found in model frame"))
      effex <- mf[, offnames, drop=FALSE]
      effect <- effect + apply(effex, 1, sum)
      #
      # construct the corresponding function
      gd <- getglmdata(model)
      goodrow <- min(which(complete.cases(gd)))
      defaultdata <- gd[goodrow, , drop=FALSE]
      effectfundata.off <- list(covname=covname,
                                fmla = formula(model),
                                defaultdata = defaultdata,
                                offnames = offnames)
      effectFun.off <- function(x) {
        d <- effectfundata.off
        # replicate default data to correct length
        df <- as.data.frame(lapply(d$defaultdata, rep, length(x)))
        # overwrite value of covariate with new data
        df[,covname] <- x
        # construct model FRAME
        mf <- model.frame(d$fmla, df)
        # check it conforms to expected structure
        if(!all(d$offnames %in% colnames(mf))) 
          stop("Internal error: mismatch in term names in effectFun")
        moff <- mf[, d$offnames, drop=FALSE]
        y <- apply(moff, 1, sum)
        return(y)
      }
    }
    if(length(termnames) == 0) {
      # Sanity clause
      # (everyone knows there ain't no Sanity Clause...)
      warning(paste("Internal error: could not find any",
                    "canonical covariates or offset terms",
                    "that depended on the covariate", sQuote(covname)))
      # Assume it's an external covariate (i.e. not used in fitted model)
      covtype <- "external"
      beta <- 0
      effect <- beta * covvalues
      effectFun <- function(x) { 0 * x }
      isoffset <- FALSE
      names(isoffset) <- covname
    }
  }

  #### Canonical covariates and coefficients
  switch(covtype,
         original={
           cancovs <- termnames
           canbeta <- termbetas
         },
         offset = ,
         canonical={
           cancovs <- covname
           canbeta <- beta
         },
         external={
           cancovs <- canbeta <- NA
         })
  
  #################################################################
  # Validate covariate values

  # locations that must have finite values 
  operative <- if(bw.restrict) insubregion & subQset else subQset

  nbg.cov <- !is.finite(covvalues)
  if(any(offending <- nbg.cov & operative)) {
    warning(paste(sum(offending), "out of", length(offending),
                  "covariate values discarded because",
                  ngettext(sum(offending), "it is", "they are"),
                  "NA or infinite"))
  }

  nbg.eff <- !is.finite(effect)
  if(any(offending <- nbg.eff & operative)) {
    warning(paste(sum(offending), "out of", length(offending),
                  "values of fitted effect discarded because",
                  ngettext(sum(offending), "it is", "they are"),
                  "NA or infinite"))
  }
  
  #################################################################
  # Restrict data to 'operative' points
  #                            with finite values
  
  nbg <- nbg.cov | nbg.eff
  ok <- !nbg & operative
  
  Q           <- Q[ok]
  covvalues   <- covvalues[ok]
  quadpoints  <- quadpoints[ok]
  resid       <- resid[ok]
  lam         <- lam[ok]
  effect      <- effect[ok]
  insubregion <- insubregion[ok]
  Z           <- Z[ok]
  wts         <- wts[ok]

  ####################################################
  # assemble data for smoothing 
  x <- covvalues
  y <- resid/wts
  if(smooth.effect) y <- y + effect 
  w <- wts
  #
  if(makefrom <- is.null(from))
    from <- min(x)
  if(maketo <- is.null(to))
    to   <- max(x)

  ####################################################
  # determine smoothing bandwidth
  #     from 'operative' data

  switch(bw.input,
         quad = {
           # bandwidth selection from covariate values at all quadrature points
           numer <- unnormdensity(x, weights=w*y,
                                  bw=bw, adjust=adjust,
                                  n=n,from=from,to=to, ...)
           sigma <- numer$bw
         },
         points= {
           # bandwidth selection from covariate values at data points
           fake <- unnormdensity(x[Z], weights=1/lam[Z],
                                 bw=bw, adjust=adjust,
                                 n=n,from=from,to=to, ...)
           sigma <- fake$bw
           numer <- unnormdensity(x, weights=w*y,
                                  bw=sigma, adjust=1,
                                  n=n,from=from,to=to, ...)
         })


  ####################################################
  # Restrict data and recompute numerator if required

  if(!is.null(subregion) && !bw.restrict) {
    # Bandwidth was computed on all data
    # Restrict to subregion and recompute numerator
    x   <- x[insubregion]
    y   <- y[insubregion]
    w   <- w[insubregion]
    Z   <- Z[insubregion]
    lam <- lam[insubregion]
    if(makefrom) from <- min(x)
    if(maketo)     to <- max(x)
    numer <- unnormdensity(x, weights=w*y,
                           bw=sigma, adjust=1,
                           n=n,from=from,to=to, ...)
  }

  ####################################################
  # Compute denominator

  denom <- unnormdensity(x, weights=w,
                         bw=sigma, adjust=1,
                         n=n,from=from,to=to, ...)

  
  ####################################################
  # Determine recommended plot range

  xr <- range(x[Z], finite=TRUE)
  alim <- xr + 0.1 * diff(xr) * c(-1,1)
  alim <- intersect.ranges(alim, c(from, to))
  
  ####################################################
  # Compute terms 

  interpolate <- function(x,y) {
    if(inherits(x, "density") && missing(y))
      approxfun(x$x, x$y, rule=2)
    else 
      approxfun(x, y, rule=2)
  }
  numfun <- interpolate(numer)
  denfun <- interpolate(denom)
  xxx <- numer$x
  yyy <- numfun(xxx)/denfun(xxx)
  # variance estimation
  # smooth 1/lambda(u) with smaller bandwidth
  tau   <- sigma/sqrt(2)
  varnumer <- unnormdensity(x, weights=w/lam,
                            bw=tau, adjust=1,
                            n=n,from=from,to=to, ...)
  varnumfun <- interpolate(varnumer)
  varestxxx <- varnumfun(xxx)/(2 * sigma * sqrt(pi) * denfun(xxx)^2)
  sd <- sqrt(varestxxx)
  # alternative estimate of variance using data points only
  varXnumer <- unnormdensity(x[Z], weights=1/lam[Z]^2,
                             bw=tau, adjust=1,
                             n=n,from=from,to=to, ...)
  varXnumfun <- interpolate(varXnumer)
  varXestxxx <- varXnumfun(xxx)/(2 * sigma * sqrt(pi) * denfun(xxx)^2)
  sdX <- sqrt(varXestxxx)
  # fitted effect
  effxxx <- effectFun(xxx)
  
  # add fitted effect of covariate, if not added before smoothing
  if(!smooth.effect)
    yyy <- yyy + effxxx
  
  ####################################################
  # pack into fv object
  
  df <- data.frame(xxx=xxx,
                   h  =yyy,
                   varh=varestxxx,
                   hi=yyy+2*sd,
                   lo=yyy-2*sd,
                   hiX=yyy+2*sdX,
                   loX=yyy-2*sdX,
                   fit=effxxx)
  # remove any funny characters in name of covariate (e.g. if it is an offset)
  Covname <- make.names(covname)
  names(df)[1] <- Covname
  desc <- c(paste("covariate", sQuote(covname)),
            "Smoothed partial residual",
            "Variance",
            "Upper limit of pointwise 5%% significance band (integral)",
            "Lower limit of pointwise 5%% significance band (integral)",
            "Upper limit of pointwise 5%% significance band (sum)",
            "Lower limit of pointwise 5%% significance band (sum)",
            paste("Parametric fitted effect of", sQuote(covname)))
  rslt <- fv(df,
             argu=Covname,
             ylab=substitute(h(X), list(X=as.name(covname))),
             valu="h",
             fmla= as.formula(paste(". ~ ", Covname)),
             alim=alim,
             labl=c(covname,
               paste("%s", paren(covname), sep=""),
               paste("var", paren(covname), sep=""),
               paste("hi", paren(covname), sep=""),
               paste("lo", paren(covname), sep=""),
               paste("hiX", paren(covname), sep=""),
               paste("loX", paren(covname), sep=""),
               paste("fit", paren(covname), sep="")),
             desc=desc,
             fname="h",
             yexp=as.expression(substitute(hat(h)(X), list(X=covname))))
  attr(rslt, "dotnames") <- c("h", "hi", "lo", "fit")
  fvnames(rslt, ".s") <- c("hi", "lo")
  # add special class data
  class(rslt) <- c("parres", class(rslt))
  attr(rslt, "stuff") <- list(covname       = paste(covname, collapse=""),
                              covtype       = covtype,
                              mediator      = mediator,
                              cancovs       = cancovs,
                              canbeta       = canbeta,
                              isoffset      = isoffset,
                              modelname     = modelname,
                              modelcall     = modelcall,
                              callstring    = callstring,
                              sigma         = sigma,
                              smooth.effect = smooth.effect,
                              restricted    = !is.null(subregion),
                              bw.input      = bw.input)
  return(rslt)
}

print.parres <- function(x, ...) {
  cat("Transformation diagnostic (class parres)\n")
  s <- attr(x, "stuff")
  cat(paste("for the", s$covtype, "covariate", sQuote(s$covname),
            if(s$covtype != "external") "in" else "for",
            "the fitted model",
            if(nchar(s$modelcall) < 30) "" else "\n\t",
            s$modelcall, "\n"))
  switch(s$covtype,
         original={
           cancovs <- s$cancovs
           med <- s$mediator
           isoffset <- s$isoffset
           if(is.null(isoffset)) isoffset <- rep.int(FALSE, length(cancovs))
           ncc <- length(cancovs)
           noff <- sum(isoffset)
           nother <- sum(!isoffset)
           explain <-
             paste(ngettext(ncc, "Fitted effect:", "Fitted effect: sum of"),
                   if(noff == 0) {
                     paste(paste(med, collapse=" and "),
                           ngettext(ncc, "term", "terms"),
                           commasep(dQuote(cancovs)))
                   } else {
                     paste(paste(med[med != "offset"], collapse=" and "),
                           ngettext(nother, "term", "terms"),
                           commasep(dQuote(cancovs[!isoffset])),
                           "and offset",
                           ngettext(noff, "term", "terms"),
                           commasep(dQuote(cancovs[isoffset])))
                   })
           cat(paste(explain, "\n"))
         },
         external={
           cat("Note: effect estimate not justified by delta method\n")
         },
         offset={},
         canonical={})
  # earlier versions were equivalent to restricted=FALSE
  if(identical(s$restricted, TRUE))
    cat("\t--Diagnostic computed for a subregion--\n")
  cat(paste("Call:", s$callstring, "\n"))
  cat(paste("Actual smoothing bandwidth sigma =", signif(s$sigma,5), "\n"))
  # earlier versions were equivalent to smooth.effect=TRUE
  sme <- !identical(s$smooth.effect, FALSE)
  if(sme) {
    cat("Algorithm: smooth(effect + residual)\n\n")
  } else {
    cat("Algorithm: effect + smooth(residual)\n\n")
  }
  NextMethod("print")
}

plot.parres <- function(x, ...) {
  xname <- deparse(substitute(x))
  do.call(plot.fv, resolve.defaults(list(x), list(...),
                                      list(main=xname, shade=c("hi", "lo"))))
}

