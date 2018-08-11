\name{latest.news}
\alias{latest.news}
\title{
  Print News About Latest Version of Package
}
\description{
  Prints the news documentation for the current version of \code{spatstat}
  or another specified package.
}
\usage{
latest.news(package = "spatstat", doBrowse=FALSE)
}
\arguments{
  \item{package}{
    Name of package for which the latest news should be printed.
  }
  \item{doBrowse}{
    Logical value indicating whether to display the results in a browser
    window instead of printing them.
  }
}
\details{
  By default, this function prints the news documentation about changes in the
  current installed version of the \pkg{spatstat} package.
  The function can be called simply by typing its
  name without parentheses (see the Examples).

  If \code{package} is given, then the function reads the
  news for the specified package from its \code{NEWS} file (if it has one) 
  and prints only the entries that refer to the current version
  of the package.

  To see the news for all previous versions as well as the current
  version, use the \R utility \code{\link[utils]{news}}. See the Examples.
}
\value{
  Null.
}
\author{
  \adrian
  and \rolf
}
\seealso{
  \code{\link[utils]{news}},
  \code{\link{bugfixes}}
}
\examples{
  if(interactive()) {

    # current news
    latest.news

    # all news
    news(package="spatstat")

  }
}
\keyword{documentation}
