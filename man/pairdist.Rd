\name{pairdist}
\alias{pairdist}
\title{Pairwise distances}
\description{
  Computes the matrix of distances between all pairs of `things'
  in a dataset
}
\usage{
  pairdist(X, \dots)
}
\arguments{
  \item{X}{
    Object specifying the locations of a set of `things'
    (such as a set of points or a set of line segments).
  }
  \item{\dots}{
    Further arguments depending on the method.
  }
}
\value{
  A square matrix whose \code{[i,j]} entry is the distance
  between the `things' numbered \code{i} and \code{j}.
}
\details{
  Given a dataset \code{X} and \code{Y}
  (representing either a point pattern or
  a line segment pattern)
  \code{pairdist} computes the distance between each pair of
  `things' in the dataset, and returns 
  a matrix containing these distances.

  The function \code{pairdist} is generic, with
  methods for point patterns (objects of class \code{"ppp"}),
  line segment patterns (objects of class \code{"psp"})
  and a default method. See the documentation
  for \code{\link{pairdist.ppp}}, \code{\link{pairdist.psp}}
  or \code{\link{pairdist.default}} for details.
}
\seealso{
  \code{\link{pairdist.ppp}},
  \code{\link{pairdist.psp}},
  \code{\link{pairdist.default}},
  \code{\link{crossdist}},
  \code{\link{nndist}},
  \code{\link{Kest}}
}
\author{Pavel Grabarnik
  \email{pavel.grabar@issp.serpukhov.su}
  and 
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{spatial}
\keyword{math}

