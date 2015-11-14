\name{intensity.ppx}
\alias{intensity.ppx}
\title{Intensity of a Multidimensional Space-Time Point Pattern}
\description{
  Calculates the intensity of points in a multi-dimensional point
  pattern of class \code{"ppx"} or \code{"pp3"}.
}
\usage{
  \method{intensity}{ppx}(X, \dots)
}
\arguments{
  \item{X}{Point pattern of class \code{"ppx"} or \code{"pp3"}.}
  \item{\dots}{Ignored.}
}
\value{
  A single number or a numeric vector.
}
\details{
  This is a method for the generic function \code{\link{intensity}}.
  It computes the empirical intensity of a multi-dimensional point pattern
  (object of class \code{"ppx"} including \code{"pp3"}),
  i.e. the average density of points per unit volume.

  If the point pattern is multitype, the intensities of the
  different types are computed separately.
}
\author{
  \adrian
  
  
  \rolf
  
  and \ege
  
}
\examples{
  X <- osteo$pts[[1]]
  intensity(X)
  marks(X) <- factor(sample(letters[1:3], npoints(X), replace=TRUE))
  intensity(X)
}





