#  First.R
#
#  $Revision: 1.32 $ $Date: 2013/04/25 06:37:43 $
#

.onLoad <- function(...) reset.spatstat.options()

.onAttach <- function(libname, pkgname) {
  store.versionstring.spatstat()
  v <- versionstring.spatstat()
  msg <- paste("spatstat", v,
               "\nType", sQuote("help(spatstat)"),
               "for an overview of spatstat",
               "\n    ", sQuote("latest.news()"),
               "for news on latest version",
               "\n    ", sQuote("licence.polygons()"),
               "for licence information on polygon calculations")
  packageStartupMessage(msg)
  invisible(NULL)
}

  
