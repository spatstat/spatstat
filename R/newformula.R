#'
#'     newformula.R
#'
#'    $Revision: 1.2 $ $Date: 2017/09/29 09:08:51 $
#' 
#'   Update formula and expand polynomial

newformula <- function(old, change, eold, enew) {
  old <- if(is.null(old)) ~1 else eval(old, eold)
  change <- if(is.null(change)) ~1 else eval(change, enew)
  old <- as.formula(old, env=eold)
  change <- as.formula(change, env=enew)
  if(spatstat.options("expand.polynom")) {
    old <- expand.polynom(old)
    change <- expand.polynom(change)
  }
  answer <- update.formula(old, change)
  return(answer)
}

