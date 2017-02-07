# Lurking variable plot for arbitrary covariate.
#
#
# $Revision: 1.52 $ $Date: 2017/02/07 08:12:05 $
#

lurking <- local({

  cumsumna <- function(x) { cumsum(ifelse(is.na(x), 0, x)) }

  ## main function
  lurking <- function(object, covariate, type="eem",
                      cumulative=TRUE,
                      clipwindow=default.clipwindow(object),
                      rv = NULL,
                      plot.sd=is.poisson(object), 
                      envelope=FALSE, nsim=39, nrank=1,
                      plot.it=TRUE,
                      typename,
                      covname, oldstyle=FALSE,
                      check=TRUE, ..., splineargs=list(spar=0.5),
                      verbose=TRUE) {
    cl <- match.call()
    ## default name for covariate
    if(missing(covname) || is.null(covname)) {
      covname <- if(is.name(cl$covariate)) as.character(cl$covariate) else
                 if(is.expression(cl$covariate)) cl$covariate else NULL
    }

    if(!identical(envelope, FALSE)) {
      ## compute simulation envelope
      Xsim <- NULL
      if(!identical(envelope, TRUE)) {
        ## some kind of object
        Y <- envelope
        if(is.list(Y) && all(sapply(Y, is.ppp))) {
          Xsim <- Y
          envelope <- TRUE
        } else if(inherits(Y, "envelope")) {
          Xsim <- attr(Y, "simpatterns")
          if(is.null(Xsim))
            stop("envelope does not contain simulated point patterns")
          envelope <- TRUE
        } else stop("Unrecognised format of argument: envelope")
        nXsim <- length(Xsim)
        if(missing(nsim) && (nXsim < nsim)) {
          warning(paste("Only", nXsim, "simulated patterns available"))
          nsim <- nXsim
        }
      }
    }
    
    ## validate object
    if(is.ppp(object)) {
      X <- object
      object <- ppm(X ~1, forcefit=TRUE)
      dont.complain.about(X)
    } else verifyclass(object, "ppm")

    ## may need to refit the model
    if(plot.sd && is.null(getglmfit(object)))
      object <- update(object, forcefit=TRUE, use.internal=TRUE)

    ## match type argument
    type <- pickoption("type", type,
                       c(eem="eem",
                         raw="raw",
                         inverse="inverse",
                         pearson="pearson",
                         Pearson="pearson"))
    if(missing(typename))
      typename <- switch(type,
                         eem="exponential energy weights",
                         raw="raw residuals",
                         inverse="inverse-lambda residuals",
                         pearson="Pearson residuals")

    ## extract spatial locations
    Q <- quad.ppm(object)
    datapoints <- Q$data
    quadpoints <- union.quad(Q)
    Z <- is.data(Q)
    wts <- w.quad(Q)
    ## subset of quadrature points used to fit model
    subQset <- getglmsubset(object)
    if(is.null(subQset)) subQset <- rep.int(TRUE, n.quad(Q))
  
    #################################################################
    ## compute the covariate
    
    if(is.im(covariate)) {
      covvalues <- covariate[quadpoints, drop=FALSE]
    } else if(is.vector(covariate) && is.numeric(covariate)) {
      covvalues <- covariate
      if(length(covvalues) != quadpoints$n)
        stop("Length of covariate vector,", length(covvalues), "!=",
             quadpoints$n, ", number of quadrature points")
    } else if(is.expression(covariate)) {
      ## Expression involving covariates in the model
      glmdata <- getglmdata(object)
      ## Fix special cases
      if(is.null(glmdata)) {
        ## default 
        glmdata <- data.frame(x=quadpoints$x, y=quadpoints$y)
        if(is.marked(quadpoints))
          glmdata$marks <- marks(quadpoints)
      }
      ## ensure x and y are in data frame 
      if(!all(c("x","y") %in% names(glmdata))) {
        glmdata$x <- quadpoints$x
        glmdata$y <- quadpoints$y
      } 
      if(!is.null(object$covariates)) {
        ## Expression may involve an external covariate that's not used in model
        neednames <- all.vars(covariate)
        if(!all(neednames %in% colnames(glmdata))) {
          moredata <- mpl.get.covariates(object$covariates, quadpoints,
                                         covfunargs=object$covfunargs)
          use <- !(names(moredata) %in% colnames(glmdata))
          glmdata <- cbind(glmdata, moredata[,use,drop=FALSE])
        }
      }
      ## Evaluate expression
      sp <- parent.frame()
      covvalues <- eval(covariate, envir= glmdata, enclos=sp)
      if(!is.numeric(covvalues))
        stop("The evaluated covariate is not numeric")
    } else 
      stop(paste("The", sQuote("covariate"), "should be either",
                 "a pixel image, an expression or a numeric vector"))

    #################################################################
    ## Validate covariate values

    nbg <- is.na(covvalues)
    if(any(offending <- nbg && subQset)) {
      if(is.im(covariate))
        warning(paste(sum(offending), "out of", length(offending),
                      "quadrature points discarded because",
                      ngettext(sum(offending), "it lies", "they lie"),
                      "outside the domain of the covariate image"))
      else
        warning(paste(sum(offending), "out of", length(offending),
                      "covariate values discarded because",
                      ngettext(sum(offending), "it is NA", "they are NA")))
    }
    ## remove points
    ok <- !nbg & subQset
    Q <- Q[ok]
    covvalues <- covvalues[ok]
    quadpoints <- quadpoints[ok]
    ## adjust
    Z <- is.data(Q)
    wts <- w.quad(Q)
    if(any(is.infinite(covvalues) | is.nan(covvalues)))
      stop("covariate contains Inf or NaN values")

    ## Quadrature points marked by covariate value
    covq <- quadpoints %mark% as.numeric(covvalues)

    ################################################################
    ## Residuals/marks attached to appropriate locations.
    ## Stoyan-Grabarnik weights are attached to the data points only.
    ## Others (residuals) are attached to all quadrature points.

    resvalues <- 
      if(!is.null(rv)) rv
      else if(type=="eem") eem(object, check=check)
      else residuals.ppm(object, type=type, check=check)
  
    if(inherits(resvalues, "msr")) {
      ## signed or vector-valued measure
      resvalues <- resvalues$val
      if(ncol(as.matrix(resvalues)) > 1)
        stop("Not implemented for vector measures; use [.msr to split into separate components")
    }

    if(type != "eem")
      resvalues <- resvalues[ok]

    res <- (if(type == "eem") datapoints else quadpoints) %mark% as.numeric(resvalues)

    ## ... and the same locations marked by the covariate
    covres <- if(type == "eem") covq[Z] else covq

    ## NAMES OF THINGS
    ## name of the covariate
    if(is.null(covname)) 
      covname <- if(is.expression(covariate)) covariate else "covariate"
    ## type of residual/mark
    if(missing(typename)) 
      typename <- if(!is.null(rv)) "rv" else ""
    
    #######################################################################
    ## START ANALYSIS
    ## Clip to subwindow if needed
    clip <-
      (!is.poisson.ppm(object) || !missing(clipwindow)) &&
      !is.null(clipwindow)
    if(clip) {
      covq <- covq[clipwindow]
      res <- res[clipwindow]
      covres <- covres[clipwindow]
      clipquad <- inside.owin(quadpoints$x, quadpoints$y, clipwindow)
      wts <- wts[ clipquad ]
    }

    ## -----------------------------------------------------------------------
    ## (A) EMPIRICAL CUMULATIVE FUNCTION
    ## based on data points if type="eem", otherwise on quadrature points

    ## Reorder the data/quad points in order of increasing covariate value
    ## and then compute the cumulative sum of their residuals/marks
    markscovres <- marks(covres)
    o <- fave.order(markscovres)
    covsort <- markscovres[o]
    cummark <- cumsumna(marks(res)[o]) 
    ## we'll plot(covsort, cummark) in the cumulative case

    ## (B) THEORETICAL MEAN CUMULATIVE FUNCTION
    ## based on all quadrature points
    
    ## Range of covariate values
    covqmarks <- marks(covq)
    covrange <- range(covqmarks, na.rm=TRUE)
    ## Suitable breakpoints
    cvalues <- seq(from=covrange[1L], to=covrange[2L], length.out=100)
    csmall <- cvalues[1L] - diff(cvalues[1:2])
    cbreaks <- c(csmall, cvalues)
    ## cumulative area as function of covariate values
    covclass <- cut(covqmarks, breaks=cbreaks)
    increm <- tapply(wts, covclass, sum)
    cumarea <- cumsumna(increm)
    ## compute theoretical mean (when model is true)
    mean0 <- if(type == "eem") cumarea else numeric(length(cumarea))
    ## we'll plot(cvalues, mean0) in the cumulative case

    ## (A'),(B') DERIVATIVES OF (A) AND (B)
    ##  Required if cumulative=FALSE  
    ##  Estimated by spline smoothing (with x values jittered)
    if(!cumulative) {
      ## fit smoothing spline to (A) 
      ss <- do.call(smooth.spline,
                    append(list(covsort, cummark),
                           splineargs)
                    )
      ## estimate derivative of (A)
      derivmark <- predict(ss, covsort, deriv=1)$y 
      ## similarly for (B) 
      ss <- do.call(smooth.spline,
                    append(list(cvalues, mean0),
                           splineargs)
                    )
      derivmean <- predict(ss, cvalues, deriv=1)$y
    }
  
    ## -----------------------------------------------------------------------
    ## Store what will be plotted
  
    if(cumulative) {
      empirical <- data.frame(covariate=covsort, value=cummark)
      theoretical <- data.frame(covariate=cvalues, mean=mean0)
    } else {
      empirical <- data.frame(covariate=covsort, value=derivmark)
      theoretical <- data.frame(covariate=cvalues, mean=derivmean)
    }

    ## ------------------------------------------------------------------------
  
    ## (C) STANDARD DEVIATION if desired
    ## (currently implemented only for Poisson)
    ## (currently implemented only for cumulative case)

    if(plot.sd && !is.poisson.ppm(object))
      warning(paste("standard deviation is calculated for Poisson model;",
                    "not valid for this model"))

    if(plot.sd && cumulative) {
      ## Fitted intensity at quadrature points
      lambda <- fitted.ppm(object, type="trend", check=check)
      lambda <- lambda[ok]
      ## Fisher information for coefficients
      asymp <- vcov(object,what="internals")
      Fisher <- asymp$fisher
      ## Local sufficient statistic at quadrature points
      suff <- asymp$suff
      suff <- suff[ok, ,drop=FALSE]
      ## Clip if required
      if(clip) {
        lambda <- lambda[clipquad]
        suff   <- suff[clipquad, , drop=FALSE]  ## suff is a matrix
      }
      ## First term: integral of lambda^(2p+1)
      switch(type,
             pearson={
               varI <- cumarea
             },
             raw={
               ## Compute sum of w*lambda for quadrature points in each interval
               dvar <- tapply(wts * lambda, covclass, sum)
               ## tapply() returns NA when the table is empty
               dvar[is.na(dvar)] <- 0
               ## Cumulate
               varI <- cumsum(dvar)
             },
             inverse=, ## same as eem
             eem={
               ## Compute sum of w/lambda for quadrature points in each interval
               dvar <- tapply(wts / lambda, covclass, sum)
               ## tapply() returns NA when the table is empty
               dvar[is.na(dvar)] <- 0
               ## Cumulate
               varI <- cumsum(dvar)
             })

      ## variance-covariance matrix of coefficients
      V <- try(solve(Fisher), silent=TRUE)
      if(inherits(V, "try-error")) {
        warning("Fisher information is singular; reverting to oldstyle=TRUE")
        oldstyle <- TRUE
      }
      
      ## Second term: B' V B
      if(oldstyle) {
        varII <- 0
      } else {
        ## lamp = lambda^(p + 1)
        lamp <- switch(type,
                       raw     = lambda, 
                       pearson = sqrt(lambda),
                       inverse =,
                       eem     = as.integer(lambda > 0))
        ## Compute sum of w * lamp * suff for quad points in intervals
        Bcontrib <- as.vector(wts * lamp) * suff
        dB <- matrix(, nrow=length(cumarea), ncol=ncol(Bcontrib))
        for(j in seq_len(ncol(dB))) 
          dB[,j] <- tapply(Bcontrib[,j], covclass, sum, na.rm=TRUE)
        ## tapply() returns NA when the table is empty
        dB[is.na(dB)] <- 0
        ## Cumulate columns
        B <- apply(dB, 2, cumsum)
        ## compute B' V B for each i 
        varII <- diag(B %*% V %*% t(B))
      }
      ##
      ## variance of residuals
      varR <- varI - varII
      ## trap numerical errors
      nbg <- (varR < 0)
      if(any(nbg)) {
        ran <- range(varR)
        varR[nbg] <- 0
        relerr <- abs(ran[1L]/ran[2L])
        nerr <- sum(nbg)
        if(relerr > 1e-6) {
          warning(paste(nerr, "negative",
                        ngettext(nerr, "value (", "values (min="),
                        signif(ran[1L], 4), ")",
                        "of residual variance reset to zero",
                        "(out of", length(varR), "values)"))
        }
      }
      theoretical$sd <- sqrt(varR)
    }

    ## 
    if(envelope) {
      ## compute envelopes by simulation
      cl$plot.it <- FALSE
      cl$envelope <- FALSE
      cl$rv <- NULL
      if(is.null(Xsim))
        Xsim <- simulate(object, nsim=nsim, progress=verbose)
      values <- NULL
      if(verbose) {
        cat("Processing.. ")
        state <- list()
      }
      for(i in seq_len(nsim)) {
        cl$object <- update(object, Xsim[[i]])
        result.i <- eval(cl, parent.frame())
        ## interpolate empirical values onto common sequence
        f.i <- with(result.i$empirical,
                    approxfun(covariate, value, rule=2))
        val.i <- f.i(theoretical$covariate)
        values <- cbind(values, val.i)
        if(verbose) state <- progressreport(i, nsim, state=state)
      }
      if(verbose) cat("Done.\n")
      hilo <- if(nrank == 1) apply(values, 1, range) else
                 apply(values, 1, orderstats, k=c(nrank, nsim-nrank+1))
      theoretical$upper <- hilo[1L,]
      theoretical$lower <- hilo[2L,]
    }
    ## ----------------  RETURN COORDINATES ----------------------------
    stuff <- list(empirical=empirical,
                  theoretical=theoretical)
    attr(stuff, "info") <- list(typename=typename,
                                cumulative=cumulative,
                                covrange=covrange,
                                covname=covname)
    class(stuff) <- "lurk"
    ## ---------------  PLOT THEM  ----------------------------------
    if(plot.it) 
      plot(stuff, ...)
    return(invisible(stuff))
  }

  lurking
})


# plot a lurk object


plot.lurk <- function(x, ..., shade="grey") {
  xplus <- append(x, attr(x, "info"))
  with(xplus, {
    ## work out plot range
    mr <- range(0, empirical$value, theoretical$mean, na.rm=TRUE)
    if(!is.null(theoretical$sd))
      mr <- range(mr,
                  theoretical$mean + 2 * theoretical$sd,
                  theoretical$mean - 2 * theoretical$sd,
                  na.rm=TRUE)
    if(!is.null(theoretical$upper))
      mr <- range(mr, theoretical$upper, theoretical$lower, na.rm=TRUE)

    ## start plot
    vname <- paste(if(cumulative)"cumulative" else "marginal", typename)
    do.call(plot,
            resolve.defaults(
              list(covrange, mr),
              list(type="n"),
              list(...),
              list(xlab=covname, ylab=vname)))
    ## Envelopes
    if(!is.null(theoretical$upper)) {
      Upper <- theoretical$upper
      Lower <- theoretical$lower
    } else if(!is.null(theoretical$sd)) {
      Upper <- with(theoretical, mean+2*sd)
      Lower <- with(theoretical, mean-2*sd)
    } else Upper <- Lower <- NULL
    if(!is.null(Upper) && !is.null(Lower)) {
      xx <- theoretical$covariate
      if(!is.null(shade)) {
        ## shaded envelope region
        shadecol <- if(is.colour(shade)) shade else "grey"
        xx <- c(xx,    rev(xx))
        yy <- c(Upper, rev(Lower))
        do.call.matched(polygon,
                        resolve.defaults(list(x=xx, y=yy),
                                         list(...),
                                         list(border=shadecol, col=shadecol)))
      } else {
        do.call(lines,
                resolve.defaults(
                  list(x = xx, y=Upper),
                  list(...),
                  list(lty=3)))
        do.call(lines,
                resolve.defaults(
                  list(x = xx, y = Lower),
                  list(...),
                  list(lty=3)))
      }
    }
    ## Empirical
    lines(value ~ covariate, empirical, ...)
    ## Theoretical mean
    do.call(lines,
            resolve.defaults(
              list(mean ~ covariate, theoretical),
              list(...),
              list(lty=2)))
  })
  return(invisible(NULL))
}



