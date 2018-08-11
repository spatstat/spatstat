#'
#'    digestCovariates.R
#'
#'     $Revision: 1.4 $  $Date: 2018/05/03 08:33:44 $
#' 

is.scov <- function(x) {
  #' Determines whether x is a valid candidate for a spatial covariate
  #' A spatial object is OK if it can be coerced to a function
  if(inherits(x, c("im", "funxy", "owin", "tess", "ssf", "leverage.ppm")))
    return(TRUE)
  #' A function(x,y,...) is OK
  if(is.function(x) && identical(names(formals(x))[1:2], c("x", "y")))
    return(TRUE)
  #' A single character "x" or "y" is OK
  if(is.character(x) && length(x) == 1 && (x %in% c("x", "y"))) 
    return(TRUE)
  #' Can't handle input
  return(FALSE)
}
  
## Assumes each input (besides W) is a single covariate or a list of covariates
## Returns a `solist` with possibly a unitname attribute
digestCovariates <- function(..., W = NULL) {
  x <- list(...)
  #' Find individual covariates in list
  valid <- sapply(x, is.scov)
  covs <- x[valid]
  #' The remaining entries are assumed to be lists of covariates
  #' so we unlist them
  x <- unlist(x[!valid], recursive = FALSE)
  valid <- sapply(x, is.scov)
  if(!all(valid))
    stop("Couldn't interpret all input as spatial covariates.")
  covs <- append(covs, x)

  if(any(needW <- !sapply(covs, is.sob))) {
    if(is.null(W)){
      boxes <- lapply(covs[!needW], Frame)
      W <- do.call(boundingbox, boxes)
    } else stopifnot(is.owin(W))
  }
  
  covunits <- vector("list", length(covs))
  # Now covs is a list of valid covariates we can loop through
  for(i in seq_along(covs)){
    covar <- covs[[i]]
    if(inherits(covar, "distfun"))
      covunits[[i]] <- unitname(covar)
    if(is.character(covar) && length(covar) == 1 && (covar %in% c("x", "y"))) {
      covar <- if(covar == "x"){
        function(x,y) { x }
      } else{
        function(x,y) { y }
      }
      covunits[[i]] <- unitname(W)
    }
    if(is.function(covar) && !inherits(covar, "funxy")){
      covar <- funxy(f = covar, W = W)
    }
    covs[[i]] <- covar
  }
  covs <- as.solist(covs)
  attr(covs, "covunits") <- covunits
  return(covs)
}
