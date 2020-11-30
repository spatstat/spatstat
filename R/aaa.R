#'
#'    aaa.R
#'
#'   Code that must be read before the rest of the R code in spatstat
#' 
#'    $Revision: 1.1 $  $Date: 2020/11/30 13:17:29 $

.spEnv <- new.env()

putSpatstatVariable <- function(name, value) {
  assign(name, value, envir=.spEnv)
}
getSpatstatVariable <- function(name) {
  get(name, envir=.spEnv)
}
existsSpatstatVariable <- function(name) {
  exists(name, envir=.spEnv)
}

