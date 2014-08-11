#
#    summary.ppm.R
#
#    summary() method for class "ppm"
#
#    $Revision: 1.65 $   $Date: 2013/11/12 16:31:54 $
#
#    summary.ppm()
#    print.summary.ppm()
#

summary.ppm <- local({
  
  covtype <- function(x) {
    if(is.im(x)) "im" else
    if(is.function(x)) "function" else
    if(is.owin(x)) "owin" else
    if(is.numeric(x) && length(x) == 1) "number" else
    if(is.factor(x)) "factor" else
    if(is.integer(x)) "integer" else
    if(is.numeric(x)) "numeric" else storage.mode(x)
  }
  
  xargs <- function(f) {
    ar <- names(formals(f))[-(1:2)]
    return(ar[ar != "..."])
  }

  summary.ppm <- function(object, ..., quick=FALSE) {
    verifyclass(object, "ppm")

    x <- object
    y <- list()
  
    #######  Extract main data components #########################

    QUAD <- object$Q
    DATA <- QUAD$data
    TREND <- x$trend

    INTERACT <- x$interaction
    if(is.null(INTERACT)) INTERACT <- Poisson()
  
    #######  Check version #########################
    
    mpl.ver <- versionstring.ppm(object)
    int.ver <- versionstring.interact(INTERACT)
    current <- versionstring.spatstat()

    virgin <- min(package_version(c(mpl.ver, int.ver)))
    
    y$antiquated <- antiquated <- (virgin <= package_version("1.5"))
    y$old        <- old        <- (virgin < majorminorversion(current))

    y$version    <- as.character(virgin)
    
    ####### Determine type of model ############################
  
    y$entries <- list()
    y$no.trend <- identical.formulae(TREND, NULL) ||
                  identical.formulae(TREND, ~1)
    y$trendvar <- trendvar <- variablesinformula(TREND)
    y$stationary <- y$no.trend || all(trendvar == "marks")

    y$poisson <- is.poisson.interact(INTERACT)

    y$marked <- is.marked.ppp(DATA)
    y$multitype <- is.multitype.ppp(DATA)
    y$marktype <- if(y$multitype) "multitype" else
                  if(y$marked) "marked" else "unmarked"

    if(y$marked) y$entries$marks <- marks(DATA)

    y$name <- paste(if(y$stationary) "Stationary " else "Nonstationary ",
                    if(y$poisson) {
                      if(y$multitype) "multitype "
                      else if(y$marked) "marked "
                      else ""
                    },
                    INTERACT$name,
                    sep="")

    ######  Fitting algorithm ########################################

    y$method <- x$method
  
    y$problems <- x$problems

    y$fitter <- if(!is.null(x$fitter)) x$fitter else "unknown"
    if(y$fitter %in% c("glm", "gam"))
      y$converged <- x$internal$glmfit$converged

    ######  Coefficients were changed after fit? #####################
  
    y$projected <- identical(x$projected, TRUE)
    y$changedcoef <- y$projected || !is.null(x$coef.orig)

    ######  Extract fitted model coefficients #########################

    y$entries$coef <- COEFS <- x$coef
    y$coef.orig <- x$coef.orig

    y$entries$Vnames <- Vnames <- x$internal$Vnames
    y$entries$IsOffset <- x$internal$IsOffset

    ###### Extract fitted interaction and summarise  #################
  
    FITIN <- fitin(x)
    y$interaction <- summary(FITIN)

    # Exit here if quick=TRUE
    
    if(identical(quick, TRUE)) {
      class(y) <- "summary.ppm"
      return(y)
    }

    ######  Does it have external covariates?  ####################

    # defaults
    y <- append(y,
                list(has.covars    = FALSE,
                     covnames      = character(0),
                     covars.used   = character(0),
                     uses.covars   = FALSE,
                     covars.are.df = FALSE,
                     expandable    = TRUE,
                     covar.type    = character(0),
                     covar.descrip = character(0),
                     has.funcs     = FALSE,
                     covfunargs    = NULL,
                     has.xargs     = FALSE,
                     xargmap       = NULL))

    if(!antiquated) {
      covars <- x$covariates
      y$has.covars <- hc <- !is.null(covars) && (length(covars) > 0)
      if(hc) {
        y$covnames <- names(covars)
        used <- (y$trendvar %in% names(covars))
        y$covars.used <- y$trendvar[used]
        y$uses.covars <- any(used)
        y$covars.are.df <- is.data.frame(covars)
        # describe covariates
        ctype <- unlist(lapply(covars, covtype))
        y$expandable <- all(ctype[used] %in%c("function", "number"))
        names(ctype) <- names(covars)
        y$covar.type <- ctype
        y$covar.descrip <- ctype
        # are there any functions?
        y$has.funcs <- any(isfun <- (ctype == "function"))
        # do covariates depend on additional arguments?
        if(y$has.funcs) {
          y$covfunargs <- x$covfunargs
          funs <- covars[isfun]
          fdescrip <- function(f) {
            if(inherits(f, "distfun")) return("distfun")
            alist <- paste(names(formals(f)), collapse=", ")
            paste("function(", alist, ")", sep="")
          }
          y$covar.descrip[isfun] <- unlist(lapply(funs, fdescrip))
          # find any extra arguments (after args 1 & 2) explicitly named
          fargs <- lapply(funs, xargs)
          nxargs <- unlist(lapply(fargs, length))
          y$has.xargs <- any(nxargs > 0)
          if(y$has.xargs) {
            # identify which function arguments are fixed in the call
            fmap <- data.frame(Covariate=rep.int(names(funs), nxargs),
                               Argument=unlist(fargs))
            fmap$Given <- (fmap$Argument %in% names(y$covfunargs))
            y$xargmap <- fmap
          }
        }
      } 
    } else {
      # Antiquated format
      # Interpret the function call instead
      callexpr <- parse(text=x$call)
      callargs <- names(as.list(callexpr[[1]]))
      # Data frame of covariates was called 'data' in versions up to 1.4-x
      y$has.covars <- !is.null(callargs) && !is.na(pmatch("data", callargs))
      # conservative guess
      y$uses.covars <- y$has.covars
      y$covfunargs <- NULL
    }
    
    ######  Arguments in call ####################################
  
    y$args <- x[c("call", "correction", "rbord")]
  
    #######  Main data components #########################

    y$entries <- append(list(quad=QUAD,
                             data=DATA,
                             interaction=INTERACT),
                        y$entries)

    if(is.character(quick) && (quick == "entries"))
      return(y)
  
    ####### Summarise data ############################

    y$data <- summary(DATA, checkdup=FALSE)
    y$quad <- summary(QUAD, checkdup=FALSE)

    if(is.character(quick) && (quick == "no prediction"))
      return(y)
  
    ######  Trend component #########################

    y$trend <- list()

    y$trend$name <- if(y$poisson) "Intensity" else "Trend"

    y$trend$formula <- if(y$no.trend) NULL else TREND

    if(y$poisson && y$no.trend) {
      # uniform Poisson process
      y$trend$value <- lambda <- exp(COEFS[[1]])
      y$trend$label <- switch(y$marktype,
                              unmarked="Uniform intensity",
                              multitype="Uniform intensity for each mark level",
                              marked="Uniform intensity in product space",
                              "")
    } else if(y$stationary) {
      # stationary
      switch(y$marktype,
             unmarked={
               # stationary non-poisson non-marked
               y$trend$label <- "First order term"
               y$trend$value <- c(beta=exp(COEFS[[1]]))
             },
             multitype={
               # stationary, multitype
               mrk <- marks(DATA)
               y$trend$label <-
                 if(y$poisson) "Intensities" else "First order terms"
               # Use predict.ppm to evaluate the fitted intensities
               lev <- factor(levels(mrk), levels=levels(mrk))
               nlev <- length(lev)
               marx <- list(x=rep.int(0, nlev), y=rep.int(0, nlev), marks=lev)
               betas <- predict(x, locations=marx, type="trend")
               names(betas) <- paste("beta_", as.character(lev), sep="")
               y$trend$value <- betas
             },
             marked={
               # stationary, marked
               y$trend$label <- "Fitted intensity coefficients"
               y$trend$value <- blankcoefnames(COEFS)
             })
    } else {
      # not stationary 
      y$trend$label <- "Fitted coefficients for trend formula"
      # extract trend terms without trying to understand them much
      if(is.null(Vnames)) 
        trendbits <- COEFS
      else {
        agree <- outer(names(COEFS), Vnames, "==")
        whichbits <- apply(!agree, 1, all)
        trendbits <- COEFS[whichbits]
      }
      y$trend$value <- blankcoefnames(trendbits)
    }
  
    # ----- parameters with SE --------------------------

    if(is.character(quick) && (quick == "no variances"))
      return(y)

    if(length(COEFS) > 0) {
      # compute standard errors
      se <- x$internal$se
      if(is.null(se)) {
        vc <- vcov(x, matrix.action="warn")
        if(!is.null(vc)) {
          se <- if(is.matrix(vc)) sqrt(diag(vc)) else
                if(length(vc) == 1) sqrt(vc) else NULL
        }
      }
      if(!is.null(se)) {
        two <- qnorm(0.975)
        lo <- COEFS - two * se
        hi <- COEFS + two * se
        pval <- 2 * pnorm(abs(COEFS)/se, lower.tail=FALSE)
        psig <- cut(pval, c(0,0.001, 0.01, 0.05, 1, Inf),
                    labels=c("***", "**", "*", "  ", "na"),
                    include.lowest=TRUE)
        notapplic <- names(COEFS) %in% c("(Intercept)", "log(lambda)")
        psig[notapplic] <- "na"
        # table of coefficient estimates with SE and 95% CI
        y$coefs.SE.CI <- data.frame(Estimate=COEFS, S.E.=se, Ztest=psig,
                                    CI95.lo=lo, CI95.hi=hi)
      }
    }
  
    class(y) <- "summary.ppm"
    return(y)
  }
  
  summary.ppm
})


coef.summary.ppm <- function(object, ...) {
  object$coefs.SE.CI
}

print.summary.ppm <- function(x, ...) {

  if(x$old)
    warning("Model was fitted by an older version of spatstat")
  
  if(is.null(x$args)) {
    # this is the quick version
    cat(paste(x$name, "\n"))
    return(invisible(NULL))
  }

  # otherwise - full details
  cat("Point process model\n")
  fitter <- if(!is.null(x$fitter)) x$fitter else "unknown"
  methodchosen <-
    if(is.null(x$method))
      "unspecified method"
    else if(fitter == "exact") "maximum likelihood" else 
      switch(x$method,
             mpl={
               if(x$poisson) {
                 # Poisson process
                 "maximum likelihood (Berman-Turner approximation)"
               } else {
                 "maximum pseudolikelihood (Berman-Turner approximation)"
               } 
             },
             logi={
               if(!x$poisson) {
                 "maximum pseudolikelihood (logistic regression approximation)"
               } else {
                 # Poisson process
                 "maximum likelihood (logistic regression approximation)"
               } 
             },
             ho="Huang-Ogata method (approximate maximum likelihood)",
             paste("unrecognised method", sQuote(x$method)))
  cat(paste("Fitting method:", methodchosen, "\n"))
  howfitted <- switch(fitter,
                      exact= "analytically",
                      gam  = "using gam()",
                      glm  = "using glm()",
                      ho   = NULL,
                      paste("using unrecognised fitter", sQuote(fitter)))
  if(!is.null(howfitted)) cat(paste("Model was fitted", howfitted, "\n"))
  if(fitter %in% c("glm", "gam")) {
    if(x$converged) cat("Algorithm converged\n")
    else cat("*** Algorithm did not converge ***\n")
  }
  if(x$projected)
    cat("Fit was projected to obtain a valid point process model\n")

  cat("Call:\n")
  print(x$args$call)

  if(x$old) 
    cat(paste("** Executed by old spatstat version", x$version, " **\n"))
  
  cat(paste("Edge correction:", dQuote(x$args$correction), "\n"))
  if(x$args$correction == "border")
    cat(paste("\t[border correction distance r =", x$args$rbord,"]\n"))

  cat("\n----------------------------------------------------\n")

  # print summary of quadrature scheme
  print(x$quad)
  
  cat("\n----------------------------------------------------\n")
  cat("FITTED MODEL:\n\n")

  # This bit is currently identical to print.ppm()
  # except for a bit more fanfare
  # and the inclusion of the 'gory details' bit
  
  notrend <-    x$no.trend
  stationary <- x$stationary
  poisson <-    x$poisson
  markeddata <- x$marked
  multitype  <- x$multitype
        
  markedpoisson <- poisson && markeddata

  # ----------- Print model type -------------------
        
  cat(x$name)
  cat("\n")

  if(markeddata) mrk <- x$entries$marks
  if(multitype) {
    cat("Possible marks: \n")
    cat(paste(levels(mrk)))
  }

  # ----- trend --------------------------

  cat(paste("\n\n ---- ", x$trend$name, ": ----\n\n", sep=""))

  if(!notrend) {
    cat("Trend formula: ")
    print(x$trend$formula)
    if(x$uses.covars) 
      cat(paste("Model depends on external",
                ngettext(length(x$covars.used), "covariate", "covariates"),
                commasep(sQuote(x$covars.used)), "\n"))
  }
  if(x$has.covars) {
    if(notrend || !x$uses.covars)
      cat("Model object contains external covariates\n")
    isdf <- identical(x$covars.are.df, TRUE)
    if(!is.null(cd <- x$covar.descrip)) {
      # print description of each covariate
      cat(paste("\nCovariates provided",
                if(isdf) " (in data frame)" else NULL,
                ":\n", sep=""))
      namescd <- names(cd)
      for(i in seq_along(cd))
        cat(paste("\t", namescd[i], ": ", cd[i], "\n", sep=""))
    }
    if(!is.null(cfa <- x$covfunargs) && length(cfa) > 0) {
      cat("Covariate function arguments (covfunargs) provided:\n")
      namescfa <- names(cfa)
      for(i in seq_along(cfa)) {
        cat(paste(namescfa[i], "= "))
        cfai <- cfa[[i]]
        if(is.numeric(cfai) && length(cfai) == 1) {
          cat(paste(cfai, "\n"))
        } else print(cfa[[i]])
      }
    }
  }

  cat(paste("\n", x$trend$label, ":\n", sep=""))
  
  tv <- x$trend$value
  if(!is.list(tv))
    print(tv)
  else 
    for(i in seq_along(tv))
      print(tv[[i]])

  # table of coefficient estimates with SE and 95% CI
  if(!is.null(cose <- x$coefs.SE.CI)) {
    cat("\n")
    print(cose)
  }
  
  # ---- Interaction ----------------------------

  if(!poisson) {
    cat("\n\n ---- Interaction: -----\n\n")
    print(x$interaction)
  }

  ####### Gory details ###################################
  cat("\n\n----------- gory details -----\n")
  COEFS <- x$entries$coef
  
  cat("\nFitted regular parameters (theta): \n")
  print(COEFS)

  cat("\nFitted exp(theta): \n")
  print(exp(unlist(COEFS)))

  ##### Warnings issued #######

  probs <- x$problems
  if(!is.null(probs) && is.list(probs) && (length(probs) > 0)) 
    lapply(probs,
           function(a) {
             if(is.list(a) && !is.null(p <- a$print))
               cat(paste("Problem:\n", p, "\n\n"))
           })
          
  return(invisible(NULL))
}

no.trend.ppm <- function(x) {
  summary.ppm(x, quick=TRUE)$no.trend
}

is.stationary <- function(x) {
  UseMethod("is.stationary")
}

is.poisson <- function(x) {
  UseMethod("is.poisson")
}

is.stationary.ppm <- function(x) {
  TREND <- x$trend
  if(is.null(TREND) || identical.formulae(TREND, ~1))
    return(TRUE)
  if(all(variablesinformula(TREND) == "marks"))
    return(TRUE)
  return(FALSE)
}

is.poisson.ppm <- function(x) {
  stopifnot(is.ppm(x))
  y <- x$interaction
  if(is.null(y)) y <- Poisson()
  is.poisson.interact(y)
}

is.marked.ppm <- function(X, ...) {
  summary.ppm(X, quick=TRUE)$marked
}

is.multitype.ppm <- function(X, ...) {
  summary.ppm(X, quick=TRUE)$multitype
}

is.expandable.ppm <- function(x) {
  return(identical(summary(x, quick="entries")$expandable, TRUE))
}

blankcoefnames <- function(x) {
  # remove name labels from ppm coefficients
  # First decide whether there are 'labels within labels'
  unlabelled <- unlist(lapply(x,
                              function(z) { is.null(names(z)) } ))
  if(all(unlabelled))
    value <- unlist(x)
  else {
    value <- list()
    for(i in seq_along(x))
      value[[i]] <- if(unlabelled[i]) unlist(x[i]) else x[[i]]
  }
  return(value)
} 
