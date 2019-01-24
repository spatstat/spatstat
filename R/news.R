#
# news.R
#
#  News and warnings
#
latest.news <- function(package="spatstat", doBrowse=FALSE, major=TRUE) {
  ## get version number
  v <- read.dcf(file=system.file("DESCRIPTION", package=package),
                fields="Version")
  if(major) {
    ## the current major version
    vp <- package_version(v)
    vv <- unlist(vp)
    v <- paste0(vv[1], ".", vv[2])
  }
  ne <- eval(substitute(news(Version >= v0, package=package), list(v0=v)))
  page(ne, method="print", doBrowse=doBrowse)
  return(invisible(ne))
}

class(latest.news) <- "autoexec"

