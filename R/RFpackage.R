#'
#'           RFpackage.R
#' 
#'    Dealing with the RandomFields package
#' 
#'    $Revision: 1.6 $  $Date: 2016/09/12 02:07:11 $

getRandomFieldsModelGen <- function(model) {
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


