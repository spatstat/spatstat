\name{requireversion}
\alias{requireversion}
\title{
  Require a Specific Version of a Package
}
\description{
  Checks that the version number of a specified package is greater than or
  equal to the specified version number.
  For use in stand-alone \R scripts.
}
\usage{
requireversion(pkg, ver, fatal=TRUE)
}
\arguments{
  \item{pkg}{
    Package name.
  }
  \item{ver}{
    Character string containing version number.
  }
  \item{fatal}{
    Logical value indicating whether an error should occur
    when the package version is less than \code{ver}.
  }
}
\details{
  This function checks whether the installed version of the
  package \code{pkg} is greater than or equal to \code{ver}.
  By default, an error occurs if this condition is not met.

  It is useful in stand-alone \R scripts, which often require
  a particular version of a package in order to work correctly.

  \bold{This function should not be used inside a package}:
  for that purpose, the dependence on packages and versions
  should be specified in the package description file.
}
\value{
  A logical value. 
}
\author{
  \adrian
}
\examples{
  requireversion(spatstat, "1.42-0")
  requireversion(spatstat, "999.999-999", fatal=FALSE)
}
\keyword{environment}
