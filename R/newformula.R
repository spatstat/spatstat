#'
#'     newformula.R
#'
#'    $Revision: 1.1 $ $Date: 2017/01/02 10:24:14 $
#' 
#'   Update formula and expand polynomial

newformula <- function(old, change, eold, enew) {
  old <- if(is.null(old)) ~1 else eval(old, eold)
  change <- if(is.null(change)) ~1 else eval(change, enew)
  old <- as.formula(old, env=eold)
  change <- as.formula(change, env=enew)
  answer <- update.formula(old, change)
  if(spatstat.options("expand.polynom")) 
    answer <- expand.polynom(answer)
  return(answer)
}

