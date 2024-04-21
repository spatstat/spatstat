#
# news.R
#
#  News and warnings
#

latest.news <- function(package=spatstat.family(), doBrowse=FALSE, major=TRUE) {
  stopifnot(is.character(package))
  ## news from spatstat.core is subsumed in spatstat.explore + spatstat.model
  defunctpackages <- "spatstat.core"
  package <- setdiff(package, defunctpackages)
  n <- length(package)
  result <- vector(mode="list", length=n)
  names(result) <- package
  for(i in seq_len(n)) {
    pack <- package[i]
    ## get version number
    v <- read.dcf(file=system.file("DESCRIPTION", package=pack),
                  fields="Version")
    if(major) {
      ## the current major version
      vp <- package_version(v)
      vv <- unlist(vp)
      v <- paste0(vv[1], ".", vv[2])
    }
    ne <- eval(substitute(news(Version >= v0, package=pack), list(v0=v)))
    if(n > 1 && !doBrowse) {
      hdr <- paste0("Package ", sQuote(pack), ":")
      if(interactive()) readline(hdr) else cat(paste(hdr, "\n"))
    }
    page(ne, method="print", doBrowse=doBrowse)
    result[[i]] <- ne
  }
  if(n == 1) result <- result[[1]]
  return(invisible(result))
}

class(latest.news) <- "autoexec"

spatstat.family <- function(subpackages=TRUE, extensions=FALSE) {
  sub <- c("spatstat.utils",
           "spatstat.data",
           "spatstat.univar",
           "spatstat.sparse",
           "spatstat.geom",
           "spatstat.random",
           "spatstat.explore",
           "spatstat.model",
           "spatstat.linnet",
           "spatstat")
  ext <- c("spatstat.gui", "spatstat.local", "spatstat.Knet")
  result <- c(if(subpackages) sub else NULL,
              if(extensions) ext else NULL)
  as.character(result)
}
