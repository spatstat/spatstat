\name{perimeter}
\Rdversion{1.1}
\alias{perimeter}
\title{
  Perimeter Length of Window
}
\description{
  Computes the perimeter length of a window
}
\usage{
perimeter(w)
}
\arguments{
  \item{w}{
    A window (object of class \code{"owin"})
    or data that can be converted to a window by \code{\link{as.owin}}.
  }
}
\details{
  This function computes the perimeter (length of the boundary) of the
  window \code{w}. If \code{w} is a rectangle or a polygonal window,
  the perimeter is the sum of the lengths of the edges of \code{w}.
  If \code{w} is a mask, it is first converted to a polygonal window
  using \code{\link{as.polygonal}}, then staircase edges are removed
  using \code{\link{simplify.owin}}, and the perimeter of the resulting
  polygon is computed.
}
\value{
  A numeric value giving the perimeter length of the window.
}
\seealso{
  \code{\link{area.owin}}
  \code{\link{diameter.owin}},
  \code{\link{owin.object}},
  \code{\link{as.owin}}
}
\examples{
  perimeter(square(3))
  data(letterR)
  perimeter(letterR)
  if(interactive()) print(perimeter(as.mask(letterR)))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
