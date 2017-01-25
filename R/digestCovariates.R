## Determines whether something is a single valid candidate for a spatial covariate
is.scov <- function(covar, warn = TRUE){
  # Spatial object is OK
  if(is.sob(covar))
    return(TRUE)
  # A function(x,y,...) is OK
  if(is.function(covar) && identical(names(formals(covar))[1:2], c("x", "y")))
    return(TRUE)
  # A single character "x" or "y" is OK
  if(is.character(covar) && length(covar) == 1 && (covar %in% c("x", "y"))) {
    return(TRUE)
  }
  # Can't handle input
  return(FALSE)
}
  
## Assumes each input (besides W) is a single covariate or a list of covariates
## Returns a `solist` with possibly a unitname attribute
digestCovariates <- function(..., W = NULL) {
  x <- list(...)
  # Find individual covariates in list
  valid <- sapply(x, is.scov)
  covs <- x[valid]
  # The remaining entries are assumed to be lists of covariates so we unlist them
  x <- unlist(x[!valid], recursive = FALSE)
  valid <- sapply(x, is.scov)
  if(!all(valid))
    stop("Couldn't interpret all input as spatial covariates.")
  covs <- append(covs, x)

  needW <- !sapply(covs, is.sob)
  if(any(needW) && is.null(W)){
    windows <- sapply(covs[!needW], as.owin, fatal = FALSE)
    W <- do.call(boundingbox, windows)
  }
  stopifnot(is.owin(W))
  
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
    if(is.function(covar)){
      covar <- funxy(f = covar, W = W)
    }
    covs[[i]] <- covar
  }
  covs <- as.solist(covs)
  attr(covs, "covunits") <- covunits
  return(covs)
}