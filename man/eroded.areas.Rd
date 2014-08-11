\name{eroded.areas}
\alias{eroded.areas}
\title{Areas of Morphological Erosions}
\description{
  Computes the areas of successive morphological erosions of a window.
}
\usage{
 eroded.areas(w, r)
}
\arguments{
  \item{w}{A window.}
  \item{r}{Numeric vector of radii at which erosions will be performed.}
}
\value{
  Numeric vector, of the same length as \code{r},
  giving the areas of the successive erosions.
}
\details{
  This function computes the areas of the erosions of the window
  \code{w} by each of the radii \code{r[i]}. 

  The morphological erosion of a set \eqn{W} by a distance \eqn{r > 0}
  is the subset 
  consisting of all points \eqn{x \in W}{x in W} such that the
  distance from \eqn{x} to the boundary of \eqn{W} is greater than
  or equal to \eqn{r}. In other words it is the result of trimming
  a margin of width \eqn{r} off the set \eqn{W}.

  The argument \code{r} should be a vector of positive numbers.
  The argument \code{w} should be a window (an object of class
  \code{"owin"}, see \code{\link{owin.object}} for details)
  or can be given in any format acceptable to \code{\link{as.owin}()}.

  Unless \code{w} is a rectangle, the computation is performed
  using a pixel raster approximation.

  To compute the eroded window itself, use \code{\link{erosion}}.
}
\seealso{
  \code{\link{owin}},
  \code{\link{as.owin}},
  \code{\link{erosion}}
}
\examples{
  w <- owin(c(0,1),c(0,1))
  a <- eroded.areas(w, seq(0.01,0.49,by=0.01))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
