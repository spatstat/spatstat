#  First.R
#
#  $Revision: 1.36 $ $Date: 2013/09/23 08:37:30 $
#

.spatstat.Nickname <- "Window Cleaner"

# previous nicknames:
# 1.32-0: "Logistical Nightmare"
# 1.33-0: "Titanic Deckchair"


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

  
