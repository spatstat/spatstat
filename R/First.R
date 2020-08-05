#  First.R
#
#  $Revision: 1.49 $ $Date: 2020/08/05 08:30:17 $
#

.onLoad <- function(...) {
  reset.spatstat.options()
  umf <- system.file("doc", "umbrella.txt", package="spatstat")
  isum <- !is.null(umf) && file.exists(umf)
  putSpatstatVariable("Spatstat.Is.Umbrella", isum)
  invisible(NULL)
}

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

  
