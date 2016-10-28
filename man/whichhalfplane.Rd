\name{whichhalfplane}
\alias{whichhalfplane}
\title{
  Test Which Side of Infinite Line a Point Falls On
}
\description{
  Given an infinite line
  and a spatial point location,
  determine which side of the line the point falls on.
}
\usage{
whichhalfplane(L, x, y = NULL)
}
\arguments{
  \item{L}{
    Object of class \code{"infline"} specifying one or more
    infinite straight lines in two dimensions.
  }
  \item{x,y}{
    Arguments acceptable to \code{\link[grDevices]{xy.coords}}
    specifying the locations of the points.
  }
}
\details{
  An infinite line \eqn{L} divides the two-dimensional plane into
  two half-planes. This function returns a matrix \code{M} of logical values
  in which \code{M[i,j] = TRUE} if the \code{j}th spatial point
  lies below or to the left of the \code{i}th line.
}
\value{
  A logical matrix.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{infline}}
}
\examples{
  L <- infline(p=runif(3), theta=runif(3, max=2*pi))
  X <- runifpoint(4)
  whichhalfplane(L, X)
}
\keyword{spatial}
\keyword{manip}
