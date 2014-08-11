#  First.R
#
#  $Revision: 1.37 $ $Date: 2013/11/13 17:18:17 $
#

.onLoad <- function(...) reset.spatstat.options()

.onAttach <- function(libname, pkgname) {
  store.versionstring.spatstat()
  ver <- versionstring.spatstat()
  ni <- read.dcf(file=system.file("DESCRIPTION", package="spatstat"),
                 fields="Nickname")
  ni <- as.character(ni)
  msg <- paste("\nspatstat", ver,
               "     ",
               paren(paste("nickname:", sQuote(ni))),
               "\nFor an introduction to spatstat, type",
               sQuote("beginner"))
  packageStartupMessage(msg)
  invisible(NULL)
}

  
