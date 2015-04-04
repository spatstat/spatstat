#
# anova.mppm.R
#
# $Revision: 1.4 $ $Date: 2015/04/04 09:05:55 $
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

  
  # Extract fit objects 
  fitz <- lapply(objex, function(x) { x$Fit$FIT })

  opt <- list(test=test, dispersion=1)

  # Finally do the appropriate ANOVA
  result <- do.call(anova, append(fitz, opt))
  
  return(result)
}

