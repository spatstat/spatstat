#
# versions.R
#
# version numbers
#
# $Revision: 1.13 $  $Date: 2019/11/29 03:38:28 $
#
#####################


# Extract version string from ppm object

versionstring.ppm <- function(object) {
  verifyclass(object, "ppm")
  v <- object$version
  if(is.null(v) || !is.list(v))
    v <- list(major=1, minor=3, release=4)
  vs <- paste(v$major, ".", v$minor, "-", v$release, sep="")
  return(vs)
}

# Extract version string from interact object

versionstring.interact <- function(object) {
  verifyclass(object, "interact")
  v <- object$version
  return(v)  # NULL before 1.11-0
}

# Get version number of current spatstat installation
# This is now saved in the spatstat cache environment
# rather than read from file every time

versionstring.spatstat <- function() {
  if(!existsSpatstatVariable("SpatstatVersion"))
    store.versionstring.spatstat()    
  getSpatstatVariable("SpatstatVersion")
}

store.versionstring.spatstat <- function() {
  vs <- read.dcf(file=system.file("DESCRIPTION", package="spatstat"),
                 fields="Version")
  vs <- as.character(vs)
  putSpatstatVariable("SpatstatVersion", vs)
}


# Extract major and minor versions only.

majorminorversion <- function(v) {
  vp <- package_version(v)
  return(package_version(paste(vp$major, vp$minor, sep=".")))
}

versioncurrency.spatstat <- function(today=Sys.Date(), checkR=TRUE) {
  ## check version currency using dates
  msg <- NULL
  if(checkR && exists("getRversion") && getRversion() >= "3.2.2") {
    ## check version of R
    rv <- R.Version()
    rdate <- with(rv, ISOdate(year, month, day))
    if(today - as.Date(rdate) > 365) {
      ## R version is really old; just warn about this
      msg <- paste(rv$version.string,
                   "is more than a year old;",
                   "we strongly recommend upgrading to the latest version")
    }
  }
  if(is.null(msg)) {
    ## check version of spatstat
    descfile <- system.file("DESCRIPTION", package="spatstat")
    packdate <- as.Date(read.dcf(file=descfile, fields="Date"))
    elapsed <- today - packdate
    if(elapsed > 75) {
      if(elapsed > 365) {
        n <- floor(elapsed/365)
        unit <- "year"
        sowhat <- "we strongly recommend upgrading to the latest version."
      } else if(elapsed > 100) {
        n <- floor(elapsed/30)
        unit <- "month"
        sowhat <- "we recommend upgrading to the latest version."
      } else {
        n <- floor(elapsed/7)
        unit <- "week"
        sowhat <- "a newer version should be available."
      }
      expired <- if(n == 1) paste("a", unit) else paste(n, paste0(unit, "s"))
      ver <- versionstring.spatstat()
      msg <- paste("spatstat version", ver,
                   "is out of date by more than",
                   paste0(expired, ";"), 
                   sowhat)
    }
  }
  return(msg)
}

# legacy function

RandomFieldsSafe <- function() { TRUE }

