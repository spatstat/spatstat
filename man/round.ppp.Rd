\name{round.ppp}
\alias{round.ppp}
\alias{round.pp3}
\alias{round.ppx}
\title{
  Apply Numerical Rounding to Spatial Coordinates
}
\description{
  Apply numerical rounding
  to the spatial coordinates of a point pattern.
}
\usage{
\method{round}{ppp}(x, digits = 0)

\method{round}{pp3}(x, digits = 0)

\method{round}{ppx}(x, digits = 0)
}
\arguments{
  \item{x}{
    A spatial point pattern in any dimension
    (object of class \code{"ppp"}, \code{"pp3"} or \code{"ppx"}).
  }
  \item{digits}{
    integer indicating the number of decimal places. 
  }
}
\details{
  These functions are methods for the generic function
  \code{\link[base]{round}}.
  They apply numerical rounding to the spatial coordinates of the
  point pattern \code{x}. 
}
\value{
  A point pattern object, of the same class as \code{x}.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{rounding}} to determine whether numbers have been
  rounded.
  
  \code{\link[base]{round}} in the Base package.
}
\examples{
  round(cells, 1)
}
\keyword{spatial}
\keyword{manip}

