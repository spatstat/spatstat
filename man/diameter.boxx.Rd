\name{diameter.boxx}
\Rdversion{1.1}
\alias{diameter.boxx}
\alias{volume.boxx}
\alias{shortside.boxx}
\alias{sidelengths.boxx}
\alias{eroded.volumes.boxx}
\title{
  Geometrical Calculations for Multi-Dimensional Box
}
\description{
  Calculates the volume, diameter, shortest side,
  side lengths, or eroded volume of a multi-dimensional box.
}
\usage{
\method{diameter}{boxx}(x)

\method{volume}{boxx}(x)

\method{shortside}{boxx}(x)

\method{sidelengths}{boxx}(x)

\method{eroded.volumes}{boxx}(x, r)
}
\arguments{
  \item{x}{
    Multi-dimensional box (object of class \code{"boxx"}).
  }
  \item{r}{
    Numeric value or vector of numeric values for which eroded volumes
    should be calculated.
  }
}
\details{
  \code{diameter.boxx}, 
  \code{volume.boxx} and \code{shortside.boxx}
  compute the diameter, volume and shortest side length of the box.
  \code{sidelengths.boxx} returns the lengths of each side of the box.
  
  \code{eroded.volumes.boxx} computes, for each entry \code{r[i]},
  the volume of the smaller box obtained by removing a slab of
  thickness \code{r[i]} from each face of the box. This smaller box is
  the subset consisting of points that lie at least \code{r[i]} units
  away from the boundary of the box.
}
\value{
  For \code{diameter.boxx}, \code{shortside.boxx} and
  \code{volume.boxx}, a single numeric value.
  For \code{sidelengths.boxx}, a numeric vector of length equal to
  the number of spatial dimensions.
  For \code{eroded.volumes.boxx}, a numeric vector of the same length
  as \code{r}.  
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{boxx}}
}
\examples{
    X <- boxx(c(0,10),c(0,10),c(0,5),c(0,2))
    diameter(X) 
    volume(X)
    shortside(X)
    sidelengths(X)
    hd <- shortside(X)/2
    eroded.volumes(X, seq(0,hd, length=10))
}
\keyword{spatial}
\keyword{math}
