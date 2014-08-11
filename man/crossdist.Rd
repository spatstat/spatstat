\name{crossdist}
\alias{crossdist}
\title{Pairwise distances}
\description{
  Computes the distances between pairs of `things'
  taken from two different datasets.
}
\usage{
  crossdist(X, Y, \dots)
}
\arguments{
  \item{X,Y}{
    Two objects of the same class.
  }
  \item{\dots}{
    Additional arguments depending on the method.
  }
}
\value{
  A matrix whose \code{[i,j]} entry is the distance
  from the \code{i}-th thing in the first dataset
  to the \code{j}-th thing in the second dataset.
}
\details{
  Given two datasets \code{X} and \code{Y}
  (representing either two point patterns or
  two line segment patterns)
  \code{crossdist} computes the Euclidean distance from each thing
  in the first dataset to each thing in the second dataset,
  and returns a matrix containing these distances.

  The function \code{crossdist} is generic, with
  methods for point patterns (objects of class \code{"ppp"}),
  line segment patterns (objects of class \code{"psp"}),
  and a default method. See the documentation for
  \code{\link{crossdist.ppp}},
  \code{\link{crossdist.psp}} or
  \code{\link{crossdist.default}} for further details.
}
\seealso{
  \code{\link{crossdist.ppp}},
  \code{\link{crossdist.psp}},
  \code{\link{crossdist.default}},
  \code{\link{pairdist}},
  \code{\link{nndist}}
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{spatial}

\keyword{math}
