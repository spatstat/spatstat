#'
#'           otherpackages.R
#' 
#'    Dealing with other packages
#' 
#'    $Revision: 1.13 $  $Date: 2016/10/12 01:28:18 $

kraeverRandomFields <- function() {
  kraever("RandomFieldsUtils")
  kraever("RandomFields")
# should no longer be needed:  
#  capture.output(RandomFieldsUtils:::.onLoad())
#  capture.output(RandomFields:::.onLoad())
  return(invisible(NULL))
}

# require a namespace and optionally check whether it is attached
kraever <- function(package, fatal=TRUE) {
  if(!requireNamespace(package, quietly=TRUE)) {
    if(fatal)
      stop(paste("The package", sQuote(package), "is required"),
           call.=FALSE)
    return(FALSE)
  }
  if(spatstat.options(paste("check", package, "loaded", sep=".")) &&
    !isNamespaceLoaded(package)){
    if(fatal)
      stop(paste("The package", sQuote(package),
                 "must be loaded: please type",
                 sQuote(paste0("library", paren(package)))),
           call.=FALSE)
    return(FALSE)
  }
  return(TRUE)
}

getRandomFieldsModelGen <- function(model) {
  kraeverRandomFields()
  if(inherits(model, "RMmodelgenerator"))
    return(model)
  if(!is.character(model))
    stop(paste("'model' should be a character string",
               "or one of the functions in the RandomFields package",
               "with a name beginning 'RM'"),
         call.=FALSE)
  switch(model,
         cauchy    = RandomFields::RMcauchy,
         exponential = ,
         exp       = RandomFields::RMexp,
         gencauchy = RandomFields::RMgencauchy,
         gauss     = RandomFields::RMgauss,
         gneiting  = RandomFields::RMgneiting,
         matern    = RandomFields::RMmatern,
         nugget    = RandomFields::RMnugget,
         spheric   = RandomFields::RMspheric,
         stable    = RandomFields::RMstable,
         whittle   = RandomFields::RMwhittle,
         {
           modgen <- try(getExportedValue("RandomFields", 
                                          paste0("RM", model)),
                         silent=TRUE)
           if(inherits(modgen, "try-error") ||
              !inherits(modgen, "RMmodelgenerator"))
             stop(paste("Model", sQuote(model), "is not recognised"))
           modgen
         })
}

