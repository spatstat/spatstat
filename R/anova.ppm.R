#
#   anova.ppm.R
#
#  $Revision: 1.8 $   $Date: 2007/04/11 10:02:09 $
#

anova.ppm <- function(object, ..., test=NULL, override=FALSE) {

  # list of models
  objex <- append(list(object), list(...))
  if(!all(unlist(lapply(objex, is.ppm))))
    stop(paste("Arguments must all be", sQuote("ppm"), "objects"))

  # non-Poisson models?
  pois <- all(unlist(lapply(objex, is.poisson.ppm)))
  if(!pois) {
    whinge <- paste("Some of the fitted models are not Poisson processes:",
                    "p-values are not supported by any theory")
    if(override)
      warning(whinge)
    else
      stop(whinge)
  }

  # all models fitted by MPL?
  mplfit <- unlist(lapply(objex, function(x) { x$method=="mpl" }))
  if(!all(mplfit)) 
    stop(paste("Not all models fitted by maximum pseudolikelihood;",
               "comparison not possible"))

  # Extract glmfit objects 
  fitz <- lapply(objex, getglmfit)

  # Any trivial models? (uniform Poisson)
  trivial <- unlist(lapply(fitz, is.null))
  
  # check whether all non-trivial models were fitted using same method
  # (all using GLM or all using GAM)
  isgam <- unlist(lapply(fitz[!trivial],
                         function(x) { inherits(x, "gam") }))
  if(any(isgam) && !all(isgam))
    warning("Some, but not all, models were fitted with use.gam=TRUE;",
            "anova may be incorrect.",
            "It is recommended to refit all models with use.gam=TRUE.")
  usegam <- any(isgam)
  
  # Force any trivial models to be refitted using GLM or GAM
  if(any(trivial)) {
    # force them to be fitted using glm
    objex[trivial] <- lapply(objex[trivial], update.ppm,
                             forcefit=TRUE, use.gam=usegam,
                             envir=parent.frame())
    fitz[trivial] <- lapply(objex[trivial], getglmfit)
  }

  # Finally do the appropriate ANOVA
  anovafun <- if(any(isgam)) "anova.gam" else "anova.glm"
  
  result <- do.call(anovafun, append(fitz, list(test=test, dispersion=1)))
  
  return(result)
}
