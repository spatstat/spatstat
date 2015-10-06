#
# anova.mppm.R
#
# $Revision: 1.5 $ $Date: 2015/08/12 07:24:07 $
#

anova.mppm <- function(object, ..., test=NULL, override=FALSE) {

  # list of models
  objex <- append(list(object), list(...))

  # Check each model is an mppm object
  if(!all(unlist(lapply(objex, function(x) {inherits(x, "mppm")}))))
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

  # All models fitted using same method?
  fitter <- unique(unlist(lapply(objex, function(x) { x$Fit$fitter })))
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
  
  # Extract fit objects 
  fitz <- lapply(objex, function(x) { x$Fit$FIT })

  opt <- list(test=test, dispersion=1)
  if(fitter == "glmmPQL") opt <- list(test=test)

  # Finally do the appropriate ANOVA
  result <- do.call(anova, append(fitz, opt))
  
  return(result)
}

