#
# covariates.R
#
# evaluate covariates
#
#   $Revision: 1.2 $  $Date: 2013/04/25 06:37:43 $
#

evalCovariate <- function(covariate, locations) {
  # evaluate covariate of any kind at specified locations
  covvalues <-
    if(is.im(covariate)) 
      safelookup(covariate, locations)
    else if(is.function(covariate)) 
      covariate(locations$x, locations$y)
    else if(is.numeric(covariate) || is.factor(covariate)) {
      if(length(covariate) == 1)
        rep.int(covariate, length(locations$x))
      else if(length(covariate) == length(locations$x))
        covariate
      else stop("Inappropriate length for covariate vector")
    } else
  stop("Covariate should be an image, a function or a factor/numeric vector")
  return(covvalues)
}

ppmCovariates <- function(model) {
  # generate list of all covariates in ppm (excluding marks)
  stopifnot(is.ppm(model))
  co <- as.list(model$covariates)
  xy <- list(x=function(x,y){x}, y=function(x,y){y})
  coplus <- append(co, xy)
  return(as.anylist(coplus))
}

findCovariate <- function(covname, scope, scopename=NULL) {
  # find the named covariate in the given ppm object or list
  if(is.ppm(scope)) {
    covlist <- ppmCovariates(scope)
    if(missing(scopename)) scopename <- "covariates in model"
  } else if(is.list(scope)) {
    covlist <- scope
  } else stop("scope should be a named list of covariates, or a ppm object")
  if(!(covname %in% names(covlist))) 
    stop(paste("covariate", dQuote(covname), "not found",
               if(!is.null(scopename)) paste("amongst", scopename) else NULL))
  covlist[[covname]]
}

