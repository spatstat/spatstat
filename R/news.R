#
# news.R
#
#  News and warnings
#
latest.news <- function(package="spatstat") {
  # get version number
  v <- read.dcf(file=system.file("DESCRIPTION", package=package),
                fields="Version")
  eval(substitute(news(Version >= v0, package=package), list(v0=v)))
}

license.polygons <- licence.polygons <- function() {
  RShowDoc("POLYGONS.txt", type="txt", package="spatstat")
}
