\name{diameter.box3} %DontDeclareMethods
\Rdversion{1.1}
\alias{diameter.box3}
\alias{volume.box3}
\alias{shortside.box3}
\alias{sidelengths.box3}
\alias{eroded.volumes.box3}
\alias{shortside}
\alias{sidelengths}
\alias{eroded.volumes}
\title{
  Geometrical Calculations for Three-Dimensional Box
}
\description{
  Calculates the volume, diameter, shortest side, side lengths, 
  or eroded volume of a three-dimensional box.
}
\usage{
\method{diameter}{box3}(x)

\method{volume}{box3}(x)

shortside(x)
sidelengths(x)
eroded.volumes(x, r)

\method{shortside}{box3}(x)

\method{sidelengths}{box3}(x)

\method{eroded.volumes}{box3}(x, r)
}
\arguments{
  \item{x}{
    Three-dimensional box (object of class \code{"box3"}).
  }
  \item{r}{
    Numeric value or vector of numeric values for which eroded volumes
    should be calculated.
  }
}
\details{
  \code{diameter.box3} computes the diameter of the box.
  \code{volume.box3} computes the volume of the box.
  \code{shortside.box3} finds the shortest of the three side lengths
  of the box.
  \code{sidelengths.box3} returns all three side lengths
  of the box.

  \code{eroded.volumes} computes, for each entry \code{r[i]},
  the volume of the smaller box obtained by removing a slab of
  thickness \code{r[i]} from each face of the box. This smaller box is
  the subset consisting of points that lie at least \code{r[i]} units
  away from the boundary of the box.
}
\value{
  For \code{diameter.box3}, \code{shortside.box3} and
  \code{volume.box3}, a single numeric
  value. For \code{sidelengths.box3}, a vector of three numbers.
  For \code{eroded.volumes}, a numeric vector of the same length
  as \code{r}.  
}
\author{\adrian
  
  
  and \rolf
  
}
\seealso{
  \code{\link{as.box3}}
}
\examples{
    X <- box3(c(0,10),c(0,10),c(0,5))
    diameter(X) 
    volume(X)
    sidelengths(X)
    shortside(X)
    hd <- shortside(X)/2
    eroded.volumes(X, seq(0,hd, length=10))
}
\keyword{spatial}
\keyword{math}
