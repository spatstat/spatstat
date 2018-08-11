#
# versions.R
#
# version numbers
#
# $Revision: 1.11 $  $Date: 2016/02/09 04:41:31 $
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

# legacy function

RandomFieldsSafe <- function() { TRUE }
