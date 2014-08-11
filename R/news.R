#
# news.R
#
#  News and warnings
#
latest.news <- function(package="spatstat") {
  # get version number
  v <- read.dcf(file=system.file("DESCRIPTION", package=package),
                fields="Version")
  ne <- eval(substitute(news(Version >= v0, package=package), list(v0=v)))
  page(ne, method="print")
  return(invisible(ne))
}

class(latest.news) <- "autoexec"

