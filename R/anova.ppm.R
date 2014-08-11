#
#   anova.ppm.R
#
#  $Revision: 1.9 $   $Date: 2013/08/19 09:00:49 $
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
  
  # force all non-trivial models to be fitted using same method
  # (all using GLM or all using GAM)
  isgam <- unlist(lapply(fitz, function(x) { inherits(x, "gam") }))
  isglm <- unlist(lapply(fitz, function(x) { inherits(x, "glm") }))
  usegam <- any(isgam)
  if(usegam && any(isglm)) {
    warning("Some, but not all, models were fitted with use.gam=TRUE;",
            "refitting all models with use.gam=TRUE.")
    objex[isglm] <- lapply(objex[isglm], update.ppm,
                           forcefit=TRUE, use.gam=TRUE,
                           envir=parent.frame())
    fitz[isglm] <- lapply(objex[isglm], getglmfit)   
  }
  
  # Force any trivial models to be refitted using GLM or GAM
  if(any(trivial)) {
    # force them to be fitted using glm
    objex[trivial] <- lapply(objex[trivial], update.ppm,
                             forcefit=TRUE, use.gam=usegam,
                             envir=parent.frame())
    fitz[trivial] <- lapply(objex[trivial], getglmfit)
  }

  # Finally do the appropriate ANOVA
  result <- do.call("anova", append(fitz, list(test=test, dispersion=1)))
  
  return(result)
}
