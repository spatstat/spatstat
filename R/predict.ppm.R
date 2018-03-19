#
#    predict.ppm.S
#
#	$Revision: 1.104 $	$Date: 2018/03/09 02:43:15 $
#
#    predict.ppm()
#	   From fitted model obtained by ppm(),	
#	   evaluate the fitted trend or conditional intensity 
#	   at a grid/list of other locations 
#
#
# -------------------------------------------------------------------

predict.ppm <- local({
  ##
  ##  extract undocumented/outdated arguments, and trap others
  ##
  xtract <- function(..., newdata=NULL, sumobj=NULL, E=NULL, total=NULL,
                     getoutofjail=FALSE) {
    if(!is.null(newdata))
      warning(paste("The use of the argument", sQuote("newdata"),
                    "is out-of-date. See help(predict.ppm)"))
    if(!is.null(total)) 
      message(paste("The use of the argument", sQuote("total"),
                    "is out-of-date. See help(predict.ppm)"))
    trap.extra.arguments(..., .Context="In predict.ppm")
    return(list(sumobj=sumobj, E=E, total=total, getoutofjail=getoutofjail))
  }
  ##
  ## confidence/prediction intervals for number of points
  predconfPois <- function(region, object, level,
                           what=c("estimate", "se",
                             "confidence", "prediction")) {
    what <- match.arg(what)
    stopifnot(0 < level && level < 1)
    lam <- predict(object, window=region)
    mu.hat <- integral.im(lam)
    if(what == "estimate") return(mu.hat)
    mo <- model.images(object, W=as.owin(lam))
    ZL <- unlist(lapply(mo,
                        function(z, w) integral.im(eval.im(z * w)),
                        w = lam))
    ZL <- matrix(ZL, nrow=1)
    var.muhat <- as.numeric(ZL %*% vcov(object) %*% t(ZL))
    sd.muhat <- sqrt(var.muhat)
    if(what == "se") return(sd.muhat)
    alpha2 <- (1-level)/2
    pp <- sort(c(alpha2, 1-alpha2))
    out <- switch(what,
                  confidence = mu.hat + qnorm(pp) * sd.muhat,
                  prediction = qmixpois(pp, mu.hat, sd.muhat, I))
    names(out) <- paste0(signif(100 * pp, 3), "%")
    out
  }

  typepublic <- c("trend", "cif", "intensity", "count")
  typeaccept <- c(typepublic, "lambda", "se", "SE", "covariates")
  typeuse    <- c(typepublic, "cif",    "se", "se", "covariates")
  
  predict.ppm <- function(object, window=NULL, ngrid=NULL, locations=NULL,
                          covariates=NULL,
                          type=c("trend", "cif", "intensity", "count"),
                          se=FALSE,
                          interval=c("none", "confidence", "prediction"),
                          level = 0.95,
                          X=data.ppm(object),
                          correction,
                          ignore.hardcore=FALSE,
                          ...,
                          dimyx=NULL, eps=NULL, 
                          new.coef=NULL, check=TRUE, repair=TRUE) {
    interval <- match.arg(interval)
    ## extract undocumented arguments 
    xarg <- xtract(...)
    sumobj <- xarg$sumobj
    E      <- xarg$E
    total  <- xarg$total
    getoutofjail <- xarg$getoutofjail
    ## match 'type' argument including 'legacy' options
    seonly <- FALSE
    if(missing(type)) type <- type[1] else {
      if(length(type) > 1) stop("Argument 'type' should be a single value")
      mt <- pmatch(type, typeaccept)
      if(is.na(mt)) stop("Argument 'type' should be one of",
                         commasep(sQuote(typepublic), " or "))
      type <- typeuse[mt]
      if(type == "se") {
        if(!getoutofjail)
          message(paste("Outdated syntax:",
                        "type='se' should be replaced by se=TRUE;",
                        "then the standard error is predict(...)$se"))
        type <- "trend"
        se <- TRUE
        seonly <- TRUE
      }
    } 
    if(!is.null(total)) {
      message("Outdated argument 'total': use 'window' and set type='count'")
      type <- "count" 
      if(!is.logical(total))
        window <- if(is.tess(total)) total else as.owin(total)
    }
    ##
    model <- object
    verifyclass(model, "ppm")
    ##  
    if(check && damaged.ppm(object)) {
      if(!repair)
        stop("object format corrupted; try update(object, use.internal=TRUE)")
      message("object format corrupted; repairing it.")
      object <- update(object, use.internal=TRUE)
    }

    if(missing(correction) || is.null(correction))
      correction <- object$correction
  
    fitcoef <- coef(object)
    if(!is.null(new.coef)) {
      ## validate coefs
      if(length(new.coef) != length(fitcoef))
        stop(paste("Argument new.coef has wrong length",
                   length(new.coef), ": should be", length(fitcoef)))
      coeffs <- new.coef
    } else {
      coeffs <- fitcoef
    }

    ##       find out what kind of model it is
    if(is.null(sumobj))
      sumobj <- summary(model, quick="entries")  # undocumented hack!
#    stationary  <- sumobj$stationary
    poisson     <- sumobj$poisson
    marked      <- sumobj$marked
    multitype   <- sumobj$multitype
    notrend     <- sumobj$no.trend
    changedcoef <- sumobj$changedcoef || !is.null(new.coef)
    trivial     <- poisson && notrend
  
    need.covariates <- sumobj$uses.covars
    covnames.needed <- sumobj$covars.used

    if(sumobj$antiquated)
      warning("The model was fitted by an out-of-date version of spatstat")  

    ##       determine mark space
    if(marked) {
      if(!multitype)
        stop("Prediction not yet implemented for general marked point processes")
      else 
        types <- levels(marks(sumobj$entries$data))
    }

    ## For Poisson models cif=intensity=trend
    if(poisson && type %in% c("cif", "intensity"))
      type <- "trend"

    ## ............. trap un-implemented cases ...................
    
    ## Standard errors not yet available for cif, intensity
    if(se && type %in% c("cif", "intensity"))
      stop(paste("Standard error for", type, "is not yet implemented"),
           call.=FALSE)

    ## Intervals are only available for unmarked Poisson models
    if(type == "count" && interval != "none" && (marked || !poisson)) {
      stop(paste0(interval, " intervals for counts are only implemented for",
                  if(marked) " unmarked" else "",
                  if(!poisson) " Poisson",
                  " models"),
           call.=FALSE)
    }

    if(interval == "prediction" && type != "count")
      stop("Prediction intervals are only available for type='count'",
           call.=FALSE)
    
    if(interval == "confidence" && type %in% c("intensity", "cif")) 
      stop(paste("Confidence intervals are not yet available for", type),
           call.=FALSE)

    estimatename <- if(interval == "none") "estimate" else interval
    
    ## ............. start computing .............................
    
    ## Total count in a region
    
    if(type == "count") {
      ## point or interval estimate, optionally with SE
      if(is.null(window)) {
        ## domain of the original data
        if(!seonly) est <- predconfPois(NULL, model, level, estimatename)
        if(se) sem <- predconfPois(NULL, model, level, "se")
      } else if(is.tess(window)) {
        ## quadrats
        tilz <- tiles(window)
        if(!seonly) {
          est <- unlist(lapply(tilz, predconfPois,
                               object=model, level=level, what=estimatename))
          if(interval != "none") # reshape
            est <- matrix(unlist(est), byrow=TRUE, ncol=2,
                          dimnames=list(names(est), names(est[[1]])))
        }
        if(se) sem <- unlist(lapply(tilz, predconfPois,
                                    object=model, level=level, what="se"))
      } else {
        ## window
        if(!seonly) est <- predconfPois(window, model, level, estimatename)
        if(se) sem <- predconfPois(window, model, level, "se")
      }
      if(!se) return(est)
      if(seonly) return(sem)
      result <- list(est, sem)
      names(result) <- c(estimatename, "se")
      return(result)
    }

    ## .....   Predict a spatial function .......
    
    if(interval != "none") {
      ## Prepare for confidence interval 
      alpha2 <- (1-level)/2
      pp <- sort(c(alpha2, 1-alpha2))
      ci.names <- paste0(signif(100 * pp, 3), "%")
      ci.q <- qnorm(pp)
    }
    
    ##      determine what kind of output is required:
    ##      (arguments present)    (output)  
    ##         window, ngrid    ->   image
    ##         locations (mask) ->   image
    ##         locations (image) ->   image
    ##         locations (rectangle) ->  treat locations as 'window'
    ##         locations (polygonal) ->  treat locations as 'window'
    ##         locations (other) ->  data frame
    ##

    if(is.im(locations))
      locations <- as.owin(locations)
    
    if(is.null(window) && is.owin(locations) && !is.mask(locations)) {
      window <- locations
      locations <- NULL
    }

    #' incompatible:
    if(!is.null(locations)) {
      #' other arguments are incompatible
      offending <- c(!is.null(ngrid), !is.null(dimyx), !is.null(eps))
      if(any(offending)) {
        offenders <- c("grid", "dimyx", "eps")[offending]
        nbad <- sum(offending)
        stop(paste(ngettext(nbad, "The argument", "The arguments"),
                   commasep(sQuote(offenders)), 
                   ngettext(nbad, "is", "are"),
                   "incompatible with", sQuote("locations")),
             call.=FALSE)
      }
    }

    #' equivalent:
    if(!is.null(ngrid) && !is.null(dimyx))
      warning(paste("The arguments", sQuote("ngrid"), "and", sQuote("dimyx"),
                    "are equivalent: only one should be given"),
              call.=FALSE)
    
    ngrid <- ngrid %orifnull% dimyx
    
    if(is.null(ngrid) && is.null(locations)) 
      ## use regular grid
      ngrid <- rev(spatstat.options("npixel"))
    
    want.image <- is.null(locations) || is.mask(locations)
    make.grid <- !is.null(ngrid) 

    ## ##############   Determine prediction points  #####################

    if(!want.image) {
      ## (A) list of (x,y) coordinates given by `locations'
      xpredict <- locations$x
      ypredict <- locations$y
      if(is.null(xpredict) || is.null(ypredict)) {
        xy <- xy.coords(locations)
        xpredict <- xy$x
        xpredict <- xy$y
      }
      if(is.null(xpredict) || is.null(ypredict))
        stop(paste("Don't know how to extract x,y coordinates from",
                   sQuote("locations")))
      ## marks if required
      if(marked) {
        ## extract marks from data frame `locations'
        mpredict <- locations$marks 
        if(is.null(mpredict))
          stop(paste("The argument", sQuote("locations"),
                     "does not contain a column of marks",
                     "(required since the fitted model",
                     "is a marked point process)"))
        if(is.factor(mpredict)) {
          ## verify mark levels match those in model
          if(!identical(all.equal(levels(mpredict), types), TRUE)) {
            if(all(levels(mpredict) %in% types))
              mpredict <- factor(mpredict, levels=types)
            else 
              stop(paste("The marks in", sQuote("locations"),
                         "do not have the same levels as",
                         "the marks in the model"))
          }
        } else {
          ## coerce to factor if possible
          if(all(mpredict %in% types))
            mpredict <- factor(mpredict, levels=types)
          else
            stop(paste("The marks in", sQuote("locations"),
                       "do not have the same values as the marks in the model"))
        }
      }
    } else {
      ## (B) pixel grid of points
      if(!make.grid) 
        ##    (B)(i) The grid is given in `locations'
        masque <- locations
      else {
        ##    (B)(ii) We have to make the grid ourselves  
        ##    Validate ngrid
        if(!is.null(ngrid)) {
          if(!is.numeric(ngrid))
            stop("ngrid should be a numeric vector")
          ngrid <- ensure2vector(ngrid)
        }
        if(is.null(window))
          window <- sumobj$entries$data$window
        masque <- as.mask(window, dimyx=ngrid, eps=eps)
      }
      ## Hack -----------------------------------------------
      ## gam with lo() will not allow extrapolation beyond the range of x,y
      ## values actually used for the fit. Check this:
      tums <- termsinformula(model$trend)
      if(any(
             tums == "lo(x)" |
             tums == "lo(y)" |
             tums == "lo(x,y)" |
             tums == "lo(y,x)")
         ) {
        ## determine range of x,y used for fit
        gg <- model$internal$glmdata
        gxr <- range(gg$x[gg$SUBSET])
        gyr <- range(gg$y[gg$SUBSET])
        ## trim window to this range
        masque <- intersect.owin(masque, owin(gxr, gyr))
      }
      ## ------------------------------------ End Hack
      ##
      ## Finally, determine x and y vectors for grid
      rxy <- rasterxy.mask(masque, drop=TRUE)
      xpredict <- rxy$x
      ypredict <- rxy$y 
    }

    ## ################  CREATE DATA FRAME  ##########################
    ##                           ... to be passed to predict.glm()  
    ##
    ## First the x, y coordinates
  
    if(!marked) 
      newdata <- data.frame(x=xpredict, y=ypredict)
    else if(!want.image) 
      newdata <- data.frame(x=xpredict, y=ypredict, marks=mpredict)
    else {
      ## replicate
      nt <- length(types)
      np <- length(xpredict)
      xpredict <- rep.int(xpredict,nt)
      ypredict <- rep.int(ypredict,nt)
      mpredict <- rep.int(types, rep.int(np, nt))
      mpredict <- factor(mpredict, levels=types)
      newdata <- data.frame(x = xpredict,
                            y = ypredict,
                            marks=mpredict)
    }

    ## ## Next the external covariates, if any
    ##
    if(need.covariates) {
      if(is.null(covariates)) {
        ## Extract covariates from fitted model object
        ## They have to be images.
        oldcov <- model$covariates
        if(is.null(oldcov))
          stop("External covariates are required, and are not available")
        if(is.data.frame(oldcov))
          stop(paste("External covariates are required.",
                     "Prediction is not possible at new locations"))
        covariates <- oldcov
      }
      ## restrict to covariates actually required for formula
      covariates <- if(is.data.frame(covariates)) {
        covariates[,covnames.needed, drop=FALSE]
      } else covariates[covnames.needed]
      covfunargs <- model$covfunargs
      covariates.df <-
        mpl.get.covariates(covariates,
                           list(x=xpredict, y=ypredict),
                           "prediction points",
                           covfunargs)
      newdata <- cbind(newdata, covariates.df)
    }

    ## ###### Set up prediction variables ################################
    ##
    ## Provide SUBSET variable
    ##
    if(is.null(newdata$SUBSET))
      newdata$SUBSET <- rep.int(TRUE, nrow(newdata))
    ##
    ## Dig out information used in Berman-Turner device 
    ##        Vnames:     the names for the ``interaction variables''
    ##        glmdata:    the data frame used for the glm fit
    ##        glmfit:     the fitted glm object
    ##

    if(!trivial) {
      Vnames <- model$internal$Vnames
      vnameprefix <- model$internal$vnameprefix
      glmdata <- getglmdata(model)
      glmfit <- getglmfit(model)
      if(object$method=="logi")
        newdata$.logi.B <- rep(glmdata$.logi.B[1], nrow(newdata))
    }

    ## Undocumented secret exit
    if(type == "covariates")
      return(list(newdata=newdata, mask=if(want.image) masque else NULL))
             
    ## ##########  COMPUTE PREDICTION ##############################
    ##
    ##   Compute the predicted value z[i] for each row of 'newdata'
    ##   Store in a vector z and reshape it later
    ##
    ##
    ## #############################################################

    needSE <- se || (interval != "none")

    attribeauts <- list()
    
    if(trivial) {
      ## ###########  UNIFORM POISSON PROCESS #####################

      lambda <- exp(coeffs[[1]])
      if(needSE) {
        npts <- nobs(model)
        se.lambda <- lambda/sqrt(npts)
      }
      switch(interval,
             none = {
               z <- rep.int(lambda, nrow(newdata))
             },
             confidence = {
               z <- matrix(lambda + se.lambda * ci.q, 
                           byrow=TRUE,
                           nrow=nrow(newdata), ncol=2,
                           dimnames=list(NULL, ci.names))
             },
             stop("Internal error: unreached"))

      if(se) 
        zse <- rep.int(se.lambda, nrow(newdata))
    
      ## ##############################################################
    } else if((type %in% c("trend", "intensity")) || poisson) {
      ##
      ## ##########  COMPUTE TREND ###################################
      ##	
      ##   set explanatory variables to zero
      ##	
      zeroes <- numeric(nrow(newdata))    
      for(vn in Vnames)    
        newdata[[vn]] <- zeroes
      ##
      ##   predict trend
      ##
      z <- lambda <- GLMpredict(glmfit, newdata, coeffs, 
                                changecoef=changedcoef)
      ##
      if(type == "intensity") 
        z <- PoisSaddle(z, fitin(model))
      
      ##
      if(needSE) {
        ## extract variance-covariance matrix of parameters
        vc <- vcov(model)
        ## compute model matrix
        fmla <- formula(model)
#        mf <- model.frame(fmla, newdata, ..., na.action=na.pass)
#        mm <- model.matrix(fmla, mf, ..., na.action=na.pass)
        mf <- model.frame(fmla, newdata, na.action=na.pass)
        mm <- model.matrix(fmla, mf, na.action=na.pass)
        if(nrow(mm) != nrow(newdata))
          stop("Internal error: row mismatch in SE calculation")
        ## compute relative variance = diagonal of quadratic form
        vv <- quadform(mm, vc)
        ## standard error
        SE <- lambda * sqrt(vv)
        if(se) 
          zse <- SE
        if(interval == "confidence") {
          z <- lambda + outer(SE, ci.q, "*")
          colnames(z) <- ci.names
        } 
      } 
      
      ## ############################################################  
    } else if(type == "cif" || type =="lambda") {
      ## ####### COMPUTE FITTED CONDITIONAL INTENSITY ################
      ##
      ## set up arguments
      inter <- model$interaction
      if(!missing(X)) stopifnot(is.ppp(X))
      W <- as.owin(data.ppm(model))
      U <- ppp(newdata$x, y=newdata$y, window=W, check=FALSE)
      if(marked) 
        marks(U) <- newdata$marks
      ## determine which prediction points are data points
      if(is.null(E))
        E <- equalpairs(U, X, marked)
    
      ## evaluate interaction
      Vnew <- evalInteraction(X, U, E, inter, correction=correction,
                              splitInf=ignore.hardcore,
                              check=check)

      if(!ignore.hardcore) {
        ## Negative infinite values of potential signify cif = zero
        cif.equals.zero <- matrowany(Vnew == -Inf)
      } else {
        ## returned as attribute, unless vacuous
        cif.equals.zero <- attr(Vnew, "-Inf") %orifnull% logical(nrow(Vnew))
      }
      attribeauts <- c(attribeauts, list(isZero=cif.equals.zero))
    
      ## Insert the potential into the relevant column(s) of `newdata'
      if(ncol(Vnew) == 1) {
        ## Potential is real valued (Vnew is a column vector)
        ## Assign values to a column of the same name in newdata
        newdata[[Vnames]] <- as.vector(Vnew)
      ##
      } else if(is.null(avail <- colnames(Vnew))) {
        ## Potential is vector-valued (Vnew is a matrix)
        ## with unnamed components.
        ## Assign the components, in order of their appearance,
        ## to the columns of newdata labelled Vnames[1], Vnames[2],... 
        for(i in seq_along(Vnames))
          newdata[[Vnames[i] ]] <- Vnew[,i]
        ##
      } else {
        ## Potential is vector-valued (Vnew is a matrix)
        ## with named components.
        ## Match variables by name
        if(all(Vnames %in% avail)) {
          for(vn in Vnames)
            newdata[[ vn ]] <- Vnew[ , vn]
        } else if(all(Vnames %in% (Pavail <- paste0(vnameprefix, avail)))) {
          for(vn in Vnames)
            newdata[[ vn ]] <- Vnew[ , match(vn, Pavail)]
        } else
          stop(paste(
            "Internal error: unable to match names",
            "of available interaction terms",
            commasep(sQuote(avail)),
            "to required interaction terms",
            commasep(sQuote(Vnames))
            ), call.=FALSE)
      }
      ## invoke predict.glm or compute prediction
      z <- GLMpredict(glmfit, newdata, coeffs, 
                      changecoef=changedcoef)
    
      ## reset to zero if potential was zero
      if(!ignore.hardcore && any(cif.equals.zero))
        z[cif.equals.zero] <- 0
    
      ## ###############################################################    
    } else
    stop(paste("Unrecognised type", sQuote(type)))

    ## ###############################################################
    ##
    ## reshape the result
    ##
    if(!want.image) {
      if(!se) {
        z <- as.vector(z)
	attributes(z) <- c(attributes(z), attribeauts)
        out <- z
      } else if(seonly) {
        out <- as.vector(zse)
      } else {
        z <- as.vector(z)
	attributes(z) <- c(attributes(z), attribeauts)
        out <- list(z, as.vector(zse))
        names(out) <- c(estimatename, "se")
      }
    }
    else {
      ## make an image of the right shape and value
      imago <- as.im(masque, value=1.0)
      if(!marked && interval=="none") {
        ## single image
        if(!se) {
          out <- imago
          ## set entries
          out[] <- z
        } else if(seonly) {
          out <- imago
          out[] <- zse
        } else {
          est <- std <- imago
          est[] <- z
          std[] <- zse
          out <- list(est, std)
          names(out) <- c(estimatename, "se")
        }
      } else if(interval != "none") {
        ## list of 2 images for CI
        if(!seonly) {
          hi <- lo <- imago
          hi[] <- z[,1]
          lo[] <- z[,2]
          est <- solist(hi, lo)
          names(est) <- ci.names
        }
        if(se) {
          std <- imago
          std[] <- zse
        }
        if(!se) {
          out <- est
        } else if(seonly) {
          out <- std
        } else {
          out <- list(est, std)
          names(out) <- c(estimatename, "se")
        }
      } else {
        ## list of images, one for each level of marks
        out <- list()
        for(i in seq_along(types)) {
          outi <- imago
          ## set entries
          outi[] <- z[newdata$marks == types[i]]
          out[[i]] <- outi
        }
        out <- as.solist(out)
        names(out) <- as.character(types)
      }
    }
    ##  
    ##  FINISHED
    ##  
    return(out)
  }

  predict.ppm
})



####################################################################
#
# compute pointwise uncertainty of fitted intensity
#
model.se.image <- function(fit, W=as.owin(fit), ..., what="sd") {
  if(!is.poisson.ppm(fit))
    stop("Only implemented for Poisson point process models", call.=FALSE)
  what <- pickoption("option", what,
                     c(sd="sd", var="var", cv="cv", CV="cv", ce="ce", CE="ce"))
  W <- as.mask(as.owin(W))
  # variance-covariance matrix of coefficients
  vc <- vcov(fit)
  np <- dim(vc)[1]
  # extract sufficient statistic for each coefficient
  mm <- model.images(fit, W, ...)
  # compute fitted intensity 
  lam <- predict(fit, locations=W)
  # initialise resulting image
  U <- as.im(W)
  U[] <- 0
  # compute pointwise matrix product, assuming vc is symmetric
  for(i in 1:np) {
    Si <- mm[[i]]
    aii <- vc[i,i]
    U <- eval.im(U + aii * Si^2)
    if(i > 1) {
      for(j in 1:(i-1)) {
        Sj <- mm[[j]]
        aij <- vc[i,j]
        twoaij <- 2 * aij
        U <- eval.im(U + twoaij * Si * Sj)
      }
    }
  }
  # the matrix product is the relative variance (CV)
  if(what=="cv")
    return(U)
  # relative sd
  if(what=="ce") {
    U <- eval.im(sqrt(U))
    return(U)
  }
  # multiply by squared intensity to obtain variance
  U <- eval.im(U * lam^2)
  # variance
  if(what=="var")
    return(U)
  # compute SD and return
  U <- eval.im(sqrt(U))
  return(U)
}

GLMpredict <- function(fit, data, coefs, changecoef=TRUE,
                       type=c("response", "link")) {
  ok <- is.finite(coefs)
  type <- match.arg(type)
  if(!changecoef && all(ok)) {
    answer <- predict(fit, newdata=data, type=type)
  } else {
    # do it by hand
    fmla <- formula(fit)
    data$.mpl.Y <- 1
    fram <- model.frame(fmla, data=data, na.action=NULL)
    # linear predictor
    mm <- model.matrix(fmla, data=fram)
    # ensure all required coefficients are present
    coefs <- fill.coefs(coefs, colnames(mm))
    ok <- is.finite(coefs)
    #
    if(all(ok)) {
      eta <- as.vector(mm %*% coefs)
    } else {
      #' ensure 0 * anything = 0
      eta <- as.vector(mm[ , ok, drop=FALSE] %*% coefs[ok])
      for(j in which(!ok)) {
        mmj <- mm[, j]
        nonzero <- is.na(mmj) | (mmj != 0)
        if(any(nonzero))
          eta[nonzero] <- eta[nonzero] + mmj[nonzero] * coefs[j]
      }
    }
    # offset
    mo <- model.offset(fram)
    if(!is.null(mo)) {
      if(is.matrix(mo))
        mo <- apply(mo, 1, sum)
      eta <- mo + eta
    }
    switch(type,
           link = {
             answer <- eta
           },
           response = {
             linkinv <- family(fit)$linkinv
             answer <- linkinv(eta)
           })
  }
  # Convert from fitted logistic prob. to lambda for logistic fit
  if(type == "response" && family(fit)$family=="binomial")
    answer <- fit$data$.logi.B[1] * answer/(1-answer)
  return(answer)
}

# An 'equalpairs' matrix E is needed in the ppm class
# to determine which quadrature points and data points are identical
# (not just which quadrature points are data points).
# It is a two-column matrix specifying all the identical pairs.
# The first column gives the index of a data point (in the data pattern X)
# and the second column gives the corresponding index in U.

# The following function determines the equal pair information
# from the coordinates (and marks) of U and X alone;
# it should be used only if we can't figure out this information otherwise.

equalpairs <- function(U, X, marked=FALSE) {
  nn <- nncross(U, X)
  coincides <- (nn$dist == 0)
  Xind <- nn$which[coincides]
  Uind <- which(coincides)
  if(marked) {
    samemarks <- (marks(X)[Xind] == marks(U)[Uind])
    Xind <- Xind[samemarks]
    Uind <- Uind[samemarks]
  }
  return(cbind(Xind, Uind))
}

  
fill.coefs <- function(coefs, required) {
  # 'coefs' should contain all the 'required' values
  coefsname <- deparse(substitute(coefs))
  nama <- names(coefs)
  if(is.null(nama)) {
    #' names cannot be matched
    if(length(coefs) != length(required))
      stop(paste("The unnamed argument", sQuote(coefsname),
                 "has", length(coefs), "entries, but",
                 length(required), "are required"),
           call.=FALSE)
    # blithely assume they match 1-1
    names(coefs) <- required
    return(coefs)
  }
  stopifnot(is.character(required))
  if(identical(nama, required)) return(coefs)
  inject <- match(nama, required)
  if(any(notneeded <- is.na(inject))) {
    warning(paste("Internal glitch: some coefficients were not required:",
                  commasep(sQuote(nama[notneeded]))),
            call.=FALSE)
    coefs <- coefs[!notneeded]
    nama <- names(coefs)
    inject <- match(nama, required)
  }
  y <- numeric(length(required))
  names(y) <- required
  y[inject] <- coefs
  return(y)
}
 
