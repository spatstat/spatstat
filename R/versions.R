#
# versions.R
#
# version numbers
#
# $Revision: 1.6 $  $Date: 2012/04/06 10:30:42 $
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
# This is now saved in the spatstat cache environment rather than read from file every time

versionstring.spatstat <- function() {
  get("SpatstatVersion", envir = .spEnv)
}

store.versionstring.spatstat <- function() {
  vs <- read.dcf(file=system.file("DESCRIPTION", package="spatstat"),
                 fields="Version")
  vs <- as.character(vs)
  assign("SpatstatVersion", vs, envir=.spEnv)
}


# Extract major and minor versions only.

majorminorversion <- function(v) {
  vp <- package_version(v)
  return(package_version(paste(vp$major, vp$minor, sep=".")))
}

# temporary test for version of RandomFields

RandomFieldsSafe <- function() {
  # The current version of RandomFields crashes in R-devel
  if(package_version(R.Version()) <= "2.15.0") return(TRUE)
  v <- read.dcf(file=system.file("DESCRIPTION", package="RandomFields"),
                fields="Version")
  v <- package_version(as.character(v))
  # Presumably it will be fixed in the next version of RandomFields
  if(v > "2.0.54") return(TRUE)
  message("RandomFields is disabled due to a bug")
  return(FALSE)
}
