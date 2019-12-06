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
requireversion(pkg, ver)
}
\arguments{
  \item{pkg}{
    Package name.
  }
  \item{ver}{
    Character string containing version number.
  }
}
\details{
  This function checks whether the installed version of the
  package \code{pkg} is greater than or equal to \code{ver}.

  It is useful in stand-alone \R scripts, which often require
  a particular version of a package in order to work correctly.

  \bold{This function should not be used inside a package}:
  for that purpose, the dependence on packages and versions
  should be specified in the package description file.
}
\value{
  Null.
}
\author{
  \adrian
}
\examples{
  requireversion(spatstat, "1.42-0")
}
\keyword{environment}
