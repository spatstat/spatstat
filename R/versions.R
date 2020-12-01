##
##  versions.R
##
## version numbers
##
##   The string P A C K A G E N A M E  is substituted by filepp
##
## $Revision: 1.20 $  $Date: 2020/11/30 10:20:24 $
##
#####################

## see also:
##     versionstring.interact()      interact.R
##     versionstring.ppm()           ppmclass.R


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
  if(checkR) {
    ## check version of R
    if(!exists("getRversion") || getRversion() < "3.3") {
      msg <- "This version of R is very outdated; an upgrade is recommended"
    } else {
      rv <- R.Version()
      rdate <- try(
        as.Date(with(rv, ISOdate(year, month, day))),
        silent=TRUE)
      if(!inherits(rdate, "Date") || length(rdate) == 0 ||
         (today - rdate > 365)) {
        ## R version is old; just warn about this
        msg <- paste(rv$version.string %orifnull% "R version",
                     "is more than a year old;",
                     "we strongly recommend upgrading to the latest version")
      }
    } 
  }
  if(is.null(msg)) {
    ## check version of spatstat
    descfile <- system.file("DESCRIPTION", package="spatstat")
    ap <- available.packages()
    cran_version <- ap[ap[, 1] == "spatstat", 2]
    local_version <- as.character(read.dcf(file=descfile, fields="Version"))
    if (identical(cran_version, local_version)) {
        elapsed <- 0
    } else {
        packdate <- as.Date(read.dcf(file=descfile, fields="Date"))
        elapsed <- today - packdate
    }
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
      msg <- paste("spatstat", "version", ver,
                   "is out of date by more than",
                   paste0(expired, ";"), 
                   sowhat)
    }
  }
  return(msg)
}

