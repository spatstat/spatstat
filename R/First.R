##  spatstat/R/First.R

.onLoad <- function(...) {
  reset.spatstat.options()
  umf <- system.file("doc", "umbrella.txt", package="spatstat")
  isum <- !is.null(umf) && file.exists(umf)
  putSpatstatVariable("Spatstat.Is.Umbrella", isum)
  invisible(NULL)
}
.onAttach <- function(libname, pkgname) {
  vs <- read.dcf(file=system.file("DESCRIPTION", package="spatstat"),
                 fields="Version")
  vs <- as.character(vs)
  putSpatstatVariable("SpatstatVersion", vs)
  nickfile <- system.file("doc", "Nickname.txt", package="spatstat")
  ni <- scan(file=nickfile, what=character(), n=1, quiet=TRUE)
  msg <- paste("\nspatstat", vs,
               "     ",
               paren(paste("nickname:", sQuote(ni))),
               "\nFor an introduction to spatstat, type",
               sQuote("beginner"), "\n")
  packageStartupMessage(msg)
  return(invisible(NULL))
}

  
