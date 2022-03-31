\name{spatstat.family}
\alias{spatstat.family}
\title{Names of All Packages in the Spatstat Family}
\description{
  Provides the names of all the packages belonging to the \pkg{spatstat}
  family of packages.
}
\usage{
  spatstat.family(subpackages=TRUE, extensions=FALSE)
}
\arguments{
  \item{subpackages}{
    Logical value specifying whether to include sub-packages.
  }
  \item{extensions}{
    Logical value specifying whether to include extension packages.
  }
}
\value{
  Character vector of package names.
}
\details{
  This function returns a character vector containing the names
  of the packages that belong to the \pkg{spatstat} family.

  By default, only the sub-packages are listed, and not the
  extension packages.
  
  A \dQuote{sub-package} is a package which is implicitly loaded or imported
  when the command \code{library(spatstat)} is issued. Currently the
  sub-packages are:
  \itemize{
    \item \code{spatstat.utils}
    \item \code{spatstat.data}
    \item \code{spatstat.sparse}
    \item \code{spatstat.geom}
    \item \code{spatstat.random}
    \item \code{spatstat.core}
    \item \code{spatstat.linnet}
    \item \code{spatstat}
  }
  An \dQuote{extension package} is a package which must be loaded
  explicitly. The extension packages are:
  \itemize{
    \item \code{spatstat.gui}
    \item \code{spatstat.local}
    \item \code{spatstat.Knet}
  }
}
\seealso{
  \code{\link{latest.news}}
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
