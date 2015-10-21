#
# anova.mppm.R
#
# $Revision: 1.6 $ $Date: 2015/10/21 09:06:57 $
#

anova.mppm <- function(object, ..., test=NULL, override=FALSE) {

  # list of models
  objex <- append(list(object), list(...))

  # Check each model is an mppm object
  if(!all(unlist(lapply(objex, is.mppm))))
    stop(paste("Arguments must all be", sQuote("mppm"), "objects"))

  # Any non-Poisson models?
  if(!all(unlist(lapply(objex, is.poisson.mppm)))) {
    whinge <- paste("Some of the fitted models are not Poisson processes:",
                    "p-values are not supported by any theory")
    if(override)
      warning(whinge)
    else
      stop(whinge)
  }

  Fits <- lapply(objex, getElement, name="Fit")
  # All models fitted using same method?
  fitter <- unique(unlist(lapply(Fits, getElement, name="fitter")))
  if(length(fitter) > 1)
    stop(paste("Models are incompatible;",
               "they were fitted by different methods (",
               paste(fitter, collapse=", "), ")" ))

  if(fitter == "glmmPQL") {
    # anova.lme requires different format of `test' argument
    # and does not recognise 'dispersion'
    if(is.null(test))
      test <- FALSE
    else {
      stopifnot(is.character(test) && length(test) == 1)
      m<- pmatch(test, c("Chisq", "F", "Cp"))
      if(is.na(m))
        stop(paste("Unrecognised test:", test))
      if(m != 1)
        stop(paste("Test", dQuote(test),
                   "is not implemented for random effects models"))
      test <- TRUE
    }
  }
  
  # Extract glm fit objects 
  fitz <- lapply(Fits, getElement, name="FIT")

  opt <- list(test=test, dispersion=1)
  if(fitter == "glmmPQL") opt <- list(test=test)

  # Finally do the appropriate ANOVA
  result <- do.call(anova, append(fitz, opt))
  
  return(result)
}

