\name{beginner}
\alias{beginner}
\title{
  Print Introduction For Beginners
}
\description{
  Prints an introduction for beginners to the \code{spatstat} package,
  or another specified package.
}
\usage{
beginner(package = "spatstat")
}
\arguments{
  \item{package}{
    Name of package.
  }
}
\details{
  This function prints an introduction for beginners
  to the \pkg{spatstat} package.

  The function can be executed simply by typing \code{beginner}
  without parentheses.

  If the argument \code{package} is given, then the function prints the
  beginner's help file \code{BEGINNER.txt} from the specified package
  (if it has one).
}
\value{
  Null.
}
\author{\adrian
  
  
  and \rolf
  
}
\seealso{
  \code{\link{latest.news}}
}
\examples{
  beginner
}
\keyword{documentation}
