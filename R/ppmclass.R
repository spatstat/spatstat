#
#	ppmclass.R
#
#	Class 'ppm' representing fitted point process models.
#
#
#	$Revision: 2.88 $	$Date: 2013/06/17 04:01:00 $
#
#       An object of class 'ppm' contains the following:
#
#            $method           model-fitting method (currently "mpl")
#
#            $coef             vector of fitted regular parameters
#                              as given by coef(glm(....))
#
#            $trend            the trend formula
#                              or NULL 
#
#            $interaction      the interaction family 
#                              (an object of class 'interact') or NULL
#
#            $Q                the quadrature scheme used
#
#            $maxlogpl         the maximised value of log pseudolikelihood
#
#            $internal         list of internal calculation results
#
#            $correction       name of edge correction method used
#            $rbord            erosion distance for border correction (or NULL)
#
#            $the.call         the originating call to ppm()
#
#            $the.version      version of mpl() which yielded the fit
#
#
#------------------------------------------------------------------------

is.ppm <- function(x) { inherits(x, "ppm") }

print.ppm <-
function(x, ...,
         what=c("all", "model", "trend", "interaction", "se", "errors")) {

  verifyclass(x, "ppm")

  misswhat <- missing(what) 

  opts <- c("model", "trend", "interaction", "se", "errors")
  what <- match.arg(what, c("all", opts), several.ok=TRUE)
  if("all" %in% what) what <- opts

  # If SE was explicitly requested, calculate it.
  # Otherwise, do it only if the model is Poisson (by default)
  do.SE <- if(!misswhat) ("se" %in% what) else
           switch(spatstat.options("print.ppm.SE"),
                  always = TRUE,
                  never  = FALSE, 
                  poisson = {
                    is.poisson(x) &&
                    !is.null(x$fitter) && (x$fitter != "gam")
                  })
  
  s <- summary.ppm(x, quick=if(do.SE) FALSE else "no variances")
        
  notrend <-    s$no.trend
  stationary <- s$stationary
  poisson <-    s$poisson
  markeddata <- s$marked
  multitype  <- s$multitype
        
  markedpoisson <- poisson && markeddata

  # ----------- Print model type -------------------

  if("model" %in% what) {
    cat(s$name)
    cat("\n")
        
    if(markeddata) mrk <- s$entries$marks
    if(multitype) {
      cat("Possible marks: \n")
      cat(paste(levels(mrk)))
      cat("\n")
    }
  }
  
  # ----- trend --------------------------

  if("trend" %in% what) {
    
#       cat(paste("\n", s$trend$name, ":\n", sep=""))

    if(!notrend) {
      cat("\nTrend formula: ")
      print(s$trend$formula, showEnv=FALSE)
    }

    tv <- s$trend$value

    if(length(tv) == 0)
      cat("\n[No trend coefficients]\n")
    else {
      cat(paste("\n", s$trend$label, ":", sep=""))
      if(is.list(tv)) {
        cat("\n")
        for(i in seq_along(tv))
          print(tv[[i]])
      } else if(is.numeric(tv) && length(tv) == 1 && is.null(names(tv))) {
        # single number: append to end of current line
        cat("\t", paste(tv), "\n")
      } else {
        # some other format 
        cat("\n")
        print(tv)
      }
    }

    cat("\n")

    if(!is.null(cfa <- s$covfunargs) && length(cfa) > 0) {
      cat("Covariate function arguments (covfunargs) provided:\n")
      for(i in seq_along(cfa)) {
        cat(paste(names(cfa)[i], "= "))
        cfai <- cfa[[i]]
        if(is.numeric(cfai) && length(cfai) == 1) {
          cat(paste(cfai, "\n"))
        } else print(cfa[[i]])
      }
    }
  }
  
  # ---- Interaction ----------------------------

  if("interaction" %in% what) {
    if(!poisson) {
      print(s$interaction, family=FALSE)
      cat("\n")
    }
  }
  
  # ----- parameter estimates with SE and 95% CI --------------------
  if("se" %in% what) {
    if(!is.null(cose <- s$coefs.SE.CI)) {
      print(cose)
    } else if(do.SE) {
      # standard error calculation failed
      cat("Standard errors unavailable; variance-covariance matrix is singular")
    } else {
      # standard error was voluntarily omitted
      cat("For standard errors, type coef(summary(x))\n")
    }
  }
  
  # ---- Warnings issued in mpl.prepare  ---------------------

  if("errors" %in% what) {
    probs <- s$problems
    if(!is.null(probs) && is.list(probs) && (length(probs) > 0)) 
      lapply(probs,
             function(x) {
               if(is.list(x) && !is.null(p <- x$print))
                 cat(paste("Problem:\n", p, "\n\n"))
             })
    
    if(s$old)
      warning(paste("Model fitted by old spatstat version", s$version))
        
  # ---- Algorithm status ----------------------------

    fitter <- s$fitter
    converged <- s$converged
    if(!is.null(fitter) && fitter %in% c("glm", "gam") && !converged)
      cat(paste("*** Fitting algorithm for", sQuote(fitter),
                "did not converge ***\n"))
  }

  if(s$projected)
    cat("Fit was projected to obtain a valid point process model\n")
  
  return(invisible(NULL))
}

quad.ppm <- function(object, drop=FALSE) {
  if(!is.ppm(object)) {
    if(inherits(object, "kppm")) object <- object$po else
    stop("object is not of class ppm")
  }
  Q <- object$Q
  if(!drop || is.null(Q))
    return(Q)
  ok <- getglmsubset(object)
  if(is.null(ok))
    return(Q)
  return(Q[ok])
}

data.ppm <- function(object) { 
  verifyclass(object, "ppm")
  object$Q$data
}

dummy.ppm <- function(object, drop=FALSE) { 
  return(quad.ppm(object, drop=drop)$dummy)
}
  
# method for 'coef'
coef.ppm <- function(object, ...) {
  verifyclass(object, "ppm")
  object$coef
}


getglmfit <- function(object) {
  verifyclass(object, "ppm")
  glmfit <- object$internal$glmfit
  if(is.null(glmfit))
      return(NULL)
  if(object$method != "mpl")
    glmfit$coefficients <- object$coef
  return(glmfit)
}

getglmdata <- function(object, drop=FALSE) {
  verifyclass(object, "ppm")
  gd <- object$internal$glmdata
  if(!drop) return(gd)
  return(gd[getglmsubset(object), , drop=FALSE])
}

getglmsubset <- function(object) {
  gd <- object$internal$glmdata
  if(object$method=="logi")
    return(gd$.logi.ok)
  return(gd$.mpl.SUBSET)
}

getppmdatasubset <- function(object) {
  # Equivalent to getglmsubset(object)[is.data(quad.ppm(object))]
  # but also works for models fitted exactly, etc
  #
  if(object$method %in% c("mpl", "ho")) {
    sub <- getglmsubset(object)
    if(!is.null(sub)) {
      Z <- is.data(quad.ppm(object))
      return(sub[Z])
    }
  }
  X <- data.ppm(object)
  sub <- if(object$correction == "border") {
    (bdist.points(X) >= object$rbord)
  } else rep(TRUE, npoints(X))
  return(sub)
}


# ??? method for 'effects' ???

valid.ppm <- function(object) {
  verifyclass(object, "ppm")
  coeffs <- coef(object)
  # ensure all coefficients are fitted, and finite
  if(!all(is.finite(coeffs)))
    return(FALSE)
  # inspect interaction
  inte <- object$interaction
  if(is.null(inte))
    return(TRUE) # Poisson process
  # extract fitted interaction coefficients
  Vnames <- object$internal$Vnames
  IsOffset <- object$internal$IsOffset  
  Icoeffs <- coeffs[Vnames[!IsOffset]]
  # check interaction
  checker <- inte$valid
  if(is.null(checker) || !newstyle.coeff.handling(inte)) {
    warning("Internal error: unable to check validity of model")
    return(NA)
  }
  answer <- checker(Icoeffs, inte)
  return(answer)
}

project.ppm <- local({
  tracemessage <- function(depth, ...) {
    if(depth == 0) return(NULL)
    spacer <- paste(rep.int("  ", depth), collapse="")
    marker <- ngettext(depth, "trace", paste("trace", depth))
    marker <- paren(marker, "[")
    cat(paste(spacer, marker, " ", paste(...), "\n", sep=""))
  }
  leaving <- function(depth) {
    tracemessage(depth, ngettext(depth, "Returning.", "Exiting level."))
  }
  project.ppm <- function(object, ..., fatal=FALSE, trace=FALSE) {
    verifyclass(object, "ppm")
    fast <- spatstat.options("project.fast")
    # user specifies 'trace' as logical
    # but 'trace' can also be integer representing trace depth
    td <- as.integer(trace)
    trace <- (td > 0)
    tdnext <- if(trace) td+1 else 0
    if(valid.ppm(object)) {
      tracemessage(td, "Model is valid.")
      leaving(td)
      return(object)
    }
    # First ensure trend coefficients are all finite
    coeffs <- coef(object)
    # Which coefficients are trend coefficients
    coefnames  <- names(coeffs)
    internames <- object$internal$Vnames
    trendnames <- coefnames[!(coefnames %in% internames)]
    # Trend terms in trend formula
    trendterms <- attr(terms(object), "term.labels")
    # Mapping from coefficients to terms of GLM
    coef2term  <- attr(model.matrix(object), "assign")
    istrend <- (coef2term > 0) & (coefnames %in% trendnames)
    # Identify non-finite trend coefficients
    bad <- istrend & !is.finite(coeffs)
    if(!any(bad)) {
      tracemessage(td, "Trend terms are valid.")
    } else {
      nbad <- sum(bad)
      tracemessage(td,
                   "Non-finite ",
                   ngettext(nbad,
                            "coefficient for term ",
                            "coefficients for terms "),
                   commasep(sQuote(trendterms[coef2term[bad]])))
      if(fast) {
        # remove first illegal term
        firstbad <- min(which(bad))
        badterm <- trendterms[coef2term[firstbad]]
        # remove this term from model
        tracemessage(td, "Removing term ", sQuote(badterm))
        removebad <- as.formula(paste("~ . - ", badterm), env=object$callframe)
        newobject <- update(object, removebad)
        if(trace) {
          tracemessage(td, "Updated model:")
          print(newobject)
        }
        # recurse
        newobject <- project.ppm(newobject, fatal=fatal, trace=tdnext)
        # return
        leaving(td)
        return(newobject)
      } else {
        # consider all illegal terms
        bestobject <- NULL
        for(i in which(bad)) {
          badterm <- trendterms[coef2term[i]]
          # remove this term from model
          tracemessage(td, "Considering removing term ", sQuote(badterm))
          removebad <- as.formula(paste("~ . - ", badterm),
                                  env=object$callframe)
          object.i <- update(object, removebad)
          if(trace) {
            tracemessage(td, "Considering updated model:")
            print(object.i)
          }
          # recurse
          object.i <- project.ppm(object.i, fatal=fatal, trace=tdnext)
          # evaluate logPL
          logPL.i   <- logLik(object.i, warn=FALSE)
          tracemessage(td, "max log pseudolikelihood = ", logPL.i)
          # optimise
          if(is.null(bestobject) || (logLik(bestobject, warn=FALSE) < logPL.i))
            bestobject <- object.i
        }
        if(trace) {
          tracemessage(td, "Best submodel:")
          print(bestobject)
        }
        # return
        leaving(td)
        return(bestobject)
      }
    } 
    # Now handle interaction
    inte <- object$interaction
    if(is.null(inte)) {
      tracemessage(td, "No interaction to check.")
      leaving(td)
      return(object)
    }
    tracemessage(td, "Inspecting interaction terms.")
    proj <- inte$project
    if(is.null(proj)) {
      whinge <- "Internal error: interaction has no projection operator"
      if(fatal) stop(whinge) 
      warning(whinge)
      leaving(td)
      return(object)
    }
    # ensure the same edge correction is used!
    correction <- object$correction
    rbord      <- object$rbord
    # apply projection 
    coef.orig <- coeffs <- coef(object)
    Vnames   <- object$internal$Vnames
    Icoeffs  <- coeffs[Vnames]
    change <- proj(Icoeffs, inte)
    if(is.null(change)) {
      tracemessage(td, "Interaction does not need updating.")
      leaving(td)
      return(object)
    }
    tracemessage(td, "Interaction is not valid.")
    if(is.numeric(change)) {
      tracemessage(td, "Interaction coefficients updated without re-fitting.")
      # old style: 'project' returned a vector of updated coefficients
      Icoeffs <- change
      # tweak interaction coefficients
      object$coef[Vnames] <- Icoeffs
      # recompute fitted interaction
      object$fitin <- NULL
      object$fitin <- fitin(object)
    } else if(is.interact(change)) {
      # new style: 'project' returns an interaction
      if(trace) {
        tracemessage(td, "Interaction changed to:")
        print(change)
      }
      # refit the whole model 
      #      (using the same edge correction)
      #      (and the same quadrature scheme)
      newobject <- update(object, interaction=change,
                          correction=correction, rbord=rbord,
                          forcefit=TRUE,
                          envir=object$callframe)
      if(trace) {
        tracemessage(td, "Updated model:")
        print(newobject)
      }
      # recurse
      newobject <- project.ppm(newobject, fatal=fatal, trace=tdnext)
      object <- newobject
    } else if(is.list(change) && all(unlist(lapply(change, is.interact)))) {
      # new style: 'project' returns a list of candidate interactions
      nchange <- length(change)
      tracemessage(td, "Considering", nchange,
                   ngettext(nchange, "submodel", "submodels"))
      bestobject <- NULL
      for(i in seq_len(nchange)) {
        change.i <- change[[i]]
        if(trace) {
          tracemessage(td,
                       "Considering", ordinal(i), 
                       "candidate submodel, with interaction:")
          print(change.i)
        }
        # refit the whole model
        object.i <- update(object, interaction=change.i,
                           correction=correction, rbord=rbord,
                           forcefit=TRUE,
                           envir=object$callframe)
        if(trace) {
          tracemessage(td, "Considering", ordinal(i),
                       "candidate updated model:")
          print(object.i)
        }
        # recurse
        object.i <- project.ppm(object.i, fatal=fatal, trace=tdnext)
        # evaluate logPL
        logPL.i   <- logLik(object.i, warn=FALSE)
        tracemessage(td, "max log pseudolikelihood = ", logPL.i)
        # optimise
        if(is.null(bestobject) || (logLik(bestobject, warn=FALSE) < logPL.i))
          bestobject <- object.i
      }
      # end loop through submodels
      if(trace) {
        tracemessage(td, "Best submodel:")
        print(bestobject)
      }
      object <- bestobject
    } else stop("Internal error: unrecognised format of update")
    object$projected <- TRUE
    object$coef.orig  <- coef.orig
    leaving(td)
    return(object)
  }
  project.ppm
})

# more methods

logLik.ppm <- function(object, ..., warn=TRUE) {
  if(!is.poisson.ppm(object) && warn) 
    warning(paste("log likelihood is not available for non-Poisson model;",
                  "log-pseudolikelihood returned"))
  method <- object$method
  switch(method,
         mpl={
           ll <- object$maxlogpl
         },
         ho={
           # evaluate the log pseudolikelihood
           Q <- quad.ppm(object, drop=TRUE)
           Z <- is.data(Q)
           w <- w.quad(Q)
           cif <- fitted(object, type="cif", drop=TRUE)
           cifdata <- cif[Z]
           ll <- sum(log(cifdata[cifdata > 0])) - sum(w * cif)
         },
         logi={
           ll <- object$maxlogpl
         },
         stop(paste("Internal error: unrecognised ppm method:",
                    dQuote(method)))
         )
  attr(ll, "df") <- length(coef(object))
  class(ll) <- "logLik"
  return(ll)
}

formula.ppm <- function(x, ...) {
  f <- x$trend
  if(is.null(f)) f <- ~1
  return(f)
}

terms.ppm <- function(x, ...) {
  terms(formula(x), ...)
}

labels.ppm <- function(object, ...) {
  # extract fitted trend coefficients
  co <- coef(object)
  Vnames <- object$internal$Vnames
  is.trend <- !(names(co) %in% Vnames)
  # model terms
  tt <- terms(object)
  lab <- attr(tt, "term.labels")
  if(length(lab) == 0)
    return(character(0))
  # model matrix
  mm <- model.matrix(object)
  ass <- attr(mm, "assign")
  # 'ass' associates coefficients with model terms
  # except ass == 0 for the Intercept
  coef.ok <- is.finite(co)
  relevant <- (ass > 0) & is.trend
  okterms <- unique(ass[coef.ok & relevant])
  return(lab[okterms])
}

extractAIC.ppm <- function (fit, scale = 0, k = 2, ...)
{
    edf <- length(coef(fit))
    aic <- AIC(fit)
    c(edf, aic + (k - 2) * edf)
}

#
# method for model.frame

model.frame.ppm <- function(formula, ...) {
  object <- formula
  gf <- getglmfit(object)
  if(is.null(gf)) {
    warning("Model re-fitted with forcefit=TRUE")
    object <- update(object, forcefit=TRUE)
    gf <- getglmfit(object)
  }
#  gd <- getglmdata(object)
#  model.frame(gf, data=gd, ...)
  model.frame(gf, ...)
}

#
# method for model.matrix

model.matrix.ppm <- function(object, data=model.frame(object),
                             ..., keepNA=TRUE) {
  gf <- getglmfit(object)
  if(is.null(gf)) {
    warning("Model re-fitted with forcefit=TRUE")
    object <- update(object, forcefit=TRUE)
    gf <- getglmfit(object)
    if(is.null(gf))
      stop("internal error: unable to extract a glm fit")
  }
  if(!keepNA) {
    # extract model matrix of glm fit object
    # restricting to its 'subset' 
    mm <- model.matrix(gf, data, ...)
    return(mm)
  }
  # extract model matrix for all cases
  mm <- model.matrix(gf, data, ..., subset=NULL)
  cn <- colnames(mm)
  gd <- getglmdata(object, drop=FALSE)
  if(nrow(mm) != nrow(gd)) {
    # can occur if covariates include NA's or interaction is -Inf
    insubset <- getglmsubset(object)
    isna <- is.na(insubset) | !insubset
    if(sum(isna) + nrow(mm) == nrow(gd)) {
      # insert rows of NA's
      mmplus <- matrix( , nrow(gd), ncol(mm))
      mmplus[isna, ] <- NA
      mmplus[!isna, ] <- mm
      mm <- mmplus
    } else 
    stop("internal error: model matrix does not match glm data frame")
  }
  colnames(mm) <- cn
  return(mm)
}

model.images <- function(object, ...) {
  UseMethod("model.images")
}

model.images.ppm <- function(object, W=as.owin(object), ...) {
  X <- data.ppm(object)
  # make a quadscheme with a dummy point at every pixel
  Q <- pixelquad(X, W)
  # construct Berman-Turner frame
  needed <- c("trend", "interaction", "covariates", "correction", "rbord")
  bt <- do.call("bt.frame", append(list(Q), object[needed]))
  # compute model matrix
  mf <- model.frame(bt$fmla, bt$glmdata, ...)
  mm <- model.matrix(bt$fmla, mf, ...)
  # retain only the entries for dummy points (pixels)
  mm <- mm[!is.data(Q), , drop=FALSE]
  # create template image
  Z <- as.im(attr(Q, "M"))
  # make images
  imagenames <- colnames(mm)
  result <- lapply(imagenames,
                   function(nama, Z, mm) {
                     values <- mm[, nama]
                     im(values, xcol=Z$xcol, yrow=Z$yrow,
                        unitname=unitname(Z))
                   },
                   Z=Z, mm=mm)
  result <- as.listof(result)
  names(result) <- imagenames
  return(result)
}

unitname.ppm <- function(x) {
  return(unitname(x$Q))
}

"unitname<-.ppm" <- function(x, value) {
  unitname(x$Q) <- value
  return(x)
}

nobs.ppm <- function(object, ...) { npoints(data.ppm(object)) }

as.interact.ppm <- function(object) {
 verifyclass(object, "ppm")
 inte <- object$interaction
 if(is.null(inte))
   inte <- Poisson()
 return(inte)
}

as.ppm <- function(object) {
  UseMethod("as.ppm")
}

as.ppm.ppm <- function(object) {
  object
}

# method for as.owin

as.owin.ppm <- function(W, ..., from=c("points", "covariates"), fatal=TRUE) {
  if(!verifyclass(W, "ppm", fatal=fatal))
    return(NULL)
  from <- match.arg(from)
  datawin <- as.owin(data.ppm(W))
  if(from == "points")
    return(datawin)
  covs <- W$covariates
  isim <- unlist(lapply(covs, is.im))
  if(!any(isim))
    return(datawin)
  cwins <- lapply(covs[isim], as.owin)
  covwin <- do.call("intersect.owin", unname(cwins))
  result <- intersect.owin(covwin, datawin)
  return(result)
}

