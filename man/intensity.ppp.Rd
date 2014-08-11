\name{intensity.ppp}  %DontDeclareMethods
\alias{intensity.ppp}
\title{
  Empirical Intensity of Point Pattern
}
\description{
  Computes the average number of points per unit area
  in a point pattern dataset.
}
\usage{
\method{intensity}{ppp}(X, ...)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  This is a method for the generic function \code{\link{intensity}}
  It computes the empirical intensity of a point pattern
  (object of class \code{"ppp"}),
  i.e. the average density of points per unit area.

  If the point pattern is multitype, the intensities of the
  different types are computed separately.
}
\value{
  A numeric value (giving the intensity) or numeric vector
  (giving the intensity for each possible type).
}
\seealso{
  \code{\link{intensity}}, 
  \code{\link{intensity.ppm}}
}
\examples{
  intensity(cells)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{nonparametric}
