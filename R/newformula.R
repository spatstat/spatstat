#'
#'     newformula.R
#'
#'    $Revision: 1.3 $ $Date: 2020/11/30 05:01:28 $
#' 
#'   Update formula and expand polynomial

newformula <- function(old, change, eold, enew,
                       expandpoly=spatstat.options("expand.polynom")) {
  old <- if(is.null(old)) ~1 else eval(old, eold)
  change <- if(is.null(change)) ~1 else eval(change, enew)
  old <- as.formula(old, env=eold)
  change <- as.formula(change, env=enew)
  if(expandpoly) {
    old <- expand.polynom(old)
    change <- expand.polynom(change)
  }
  answer <- update.formula(old, change)
  return(answer)
}

