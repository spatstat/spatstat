\name{sidelengths.owin}
\alias{sidelengths.owin}
\alias{shortside.owin}
\title{Side Lengths of Enclosing Rectangle of a Window}
\description{
  Computes the side lengths of the (enclosing rectangle of) a window.
}
\usage{
 \method{sidelengths}{owin}(x)

 \method{shortside}{owin}(x)
}
\arguments{
  \item{x}{
    A window whose side lengths will be computed.
    Object of class \code{"owin"}.
  }
}
\value{
  For \code{sidelengths.owin}, a numeric vector of length 2
  giving the side-lengths (\eqn{x} then \eqn{y}) of the enclosing rectangle.
  For \code{shortside.owin}, a numeric value.
}
\details{
  The functions \code{shortside} and \code{sidelengths} are generic.
  The functions documented here are the methods for the class \code{"owin"}.

  \code{sidelengths.owin} computes the 
  side-lengths of the enclosing rectangle of the window \code{x}.

  For safety, both functions give a warning if the window is not a rectangle.
  To suppress the warning, first convert the window to a rectangle
  using \code{\link{as.rectangle}}.

  \code{shortside.owin} computes the minimum of the two side-lengths.
}
\seealso{
  \code{\link{shortside}},   \code{\link{sidelengths}}
  for the generic functions.
  
  \code{\link{area.owin}},
  \code{\link{diameter.owin}},
  \code{\link{perimeter}}
  for other geometric calculations on \code{"owin"} objects.
  
  \code{\link{owin}},
  \code{\link{as.owin}}.
}
\examples{
  w <- owin(c(0,2),c(-1,3))
  sidelengths(w)
  shortside(as.rectangle(letterR))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
