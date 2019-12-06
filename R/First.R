#  First.R
#
#  $Revision: 1.48 $ $Date: 2019/12/06 01:38:23 $
#

.onLoad <- function(...) reset.spatstat.options()

.onAttach <- function(libname, pkgname) {
  store.versionstring.spatstat()
  ver <- versionstring.spatstat()
  ## descfile <- system.file("DESCRIPTION", package="spatstat")
  nickfile <- system.file("doc", "Nickname.txt", package="spatstat")
  ni <- scan(file=nickfile, what=character(), n=1, quiet=TRUE)
  msg <- paste("\nspatstat", ver,
               "     ",
               paren(paste("nickname:", sQuote(ni))),
               "\nFor an introduction to spatstat, type",
               sQuote("beginner"), "\n")
  packageStartupMessage(msg)
  cur <- versioncurrency.spatstat()
  if(!is.null(cur))
    packageStartupMessage(paste("\nNote:", cur))
  invisible(NULL)
}

  
