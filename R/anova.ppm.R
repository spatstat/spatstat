#
#   anova.ppm.R
#
#  $Revision: 1.13 $   $Date: 2014/06/25 10:10:49 $
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

  expandedfrom1 <- FALSE
  if(length(objex) == 1 && inherits(object, "ippm")) {
    ## we can't rely on anova.glm to get the df right in this case
    ## so we have to re-fit explicitly
    Terms <- drop.scope(object)
    if((nT <- length(Terms)) > 0) {
      ## generate models by adding terms sequentially
      objex <- vector(mode="list", length=nT+1)
      for(n in 1:nT) {
        ## model containing terms 1, ..., n-1
        fmla <- paste(". ~ . - ", paste(Terms[n:nT], collapse=" - "))
        fmla <- as.formula(fmla)
        objex[[n]] <- update(object, fmla)
      }
      ## full model
      objex[[nT+1]] <- object
      expandedfrom1 <- TRUE
    }
  }
    
  ## all models fitted by MPL?
  mplfit <- unlist(lapply(objex, function(x) { x$method=="mpl" }))
  if(!all(mplfit)) 
    stop(paste("Not all models fitted by maximum pseudolikelihood;",
               "comparison not possible"))

  ## Extract glmfit objects 
  fitz <- lapply(objex, getglmfit)

  ## Any trivial models? (uniform Poisson)
  trivial <- unlist(lapply(fitz, is.null))
  
  ## force all non-trivial models to be fitted using same method
  ## (all using GLM or all using GAM)
  isgam <- unlist(lapply(fitz, function(x) { inherits(x, "gam") }))
  isglm <- unlist(lapply(fitz, function(x) { inherits(x, "glm") }))
  usegam <- any(isgam)
  if(usegam && any(isglm)) {
    warning("Some, but not all, models were fitted with use.gam=TRUE;",
            "refitting all models with use.gam=TRUE.")
    objex[isglm] <- lapply(objex[isglm], update.ppm,
                           forcefit=TRUE, use.gam=TRUE)
    fitz[isglm] <- lapply(objex[isglm], getglmfit)   
  }
  
  ## Force any trivial models to be refitted using GLM or GAM
  if(any(trivial)) {
    # force them to be fitted using glm
    objex[trivial] <- lapply(objex[trivial], update.ppm,
                             forcefit=TRUE, use.gam=usegam)
    fitz[trivial] <- lapply(objex[trivial], getglmfit)
  }

  ## If any models were fitted by ippm we need to correct the df
  if(any(unlist(lapply(objex, inherits, what="ippm")))) {
    nfree <- unlist(lapply(lapply(objex, logLik), attr, which="df"))
    ncanonical <- unlist(lapply(lapply(objex, coef), length))
    nextra <- nfree - ncanonical
    for(i in seq_along(fitz))
      if(nextra[i] != 0)
        fitz[[i]]$df.residual <- fitz[[i]]$df.residual - nextra[i]
  }
  
  ## Finally do the appropriate ANOVA
  result <- do.call("anova", append(fitz, list(test=test, dispersion=1)))

  ## Remove approximation-dependent columns 
  result[, "Resid. Df"] <- NULL
  result[, "Resid. Dev"] <- NULL

  ## edit header 
  if(!is.null(h <- attr(result, "heading"))) {
    ## remove .mpl.Y from formulae if present
    h <- gsub(".mpl.Y", "", h)
    ## delete GLM information if present
    h <- gsub("Model: quasi, link: log", "", h)
    h <- gsub("Response: ", "", h)
    ## remove blank lines (up to 4 consecutive blanks can occur)
    for(i in 1:5)
      h <- gsub("\n\n", "\n", h)
    ## Add explanation if we did the stepwise thing ourselves
    if(expandedfrom1)
      h <- c(h[1], "Terms added sequentially (first to last)\n", h[-1])
    ## Contract spaces in output if spatstat.options('terse') >= 2
    if(!waxlyrical('space'))
      h <- gsub("\n$", "", h)
    ## Put back
    attr(result, "heading") <- h
  }
  
  return(result)
}
