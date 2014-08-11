\name{as.data.frame.ppp}
\alias{as.data.frame.ppp}
\title{Coerce Point Pattern to a Data Frame}
\description{
  Extracts the coordinates of the points in a point pattern,
  and their marks if any, and returns them in a data frame.
}
\usage{
\method{as.data.frame}{ppp}(x, row.names = NULL, ...)
}
\arguments{
  \item{x}{Point pattern (object of class \code{"ppp"}).}
  \item{row.names}{Optional character vector of row names.}
  \item{\dots}{Ignored.}
}
\details{
  This is a method for the generic function \code{\link{as.data.frame}}
  for the class \code{"ppp"} of point patterns.

  It extracts the coordinates of the points in the point pattern, and
  returns them as columns named \code{x} and \code{y} in a data frame.
  If the points were marked, the marks are returned as a column
  named \code{marks} with the same type as in the point pattern dataset.
}
\value{
  A data frame.
}
\examples{
  data(amacrine)
  df <- as.data.frame(amacrine)
  df[1:5,]
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
