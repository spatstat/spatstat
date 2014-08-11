#  First.R
#
#  $Revision: 1.34 $ $Date: 2013/08/02 05:48:26 $
#

.spatstat.Nickname <- "Titanic Deckchair"

# previous nicknames:
# 1.32-0: "Logistical Nightmare"

.onLoad <- function(...) reset.spatstat.options()

.onAttach <- function(libname, pkgname) {
  store.versionstring.spatstat()
  ver <- versionstring.spatstat()
  nick <- .spatstat.Nickname
  msg <- paste("\nspatstat", ver,
               "     ",
               paren(paste("nickname:", sQuote(nick))),
               "\nFor an introduction to spatstat, type",
               sQuote("beginner"))
  packageStartupMessage(msg)
  invisible(NULL)
}

  
