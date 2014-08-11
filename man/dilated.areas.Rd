\name{dilated.areas}
\Rdversion{1.1}
\alias{dilated.areas}
\title{
  Areas of Morphological Dilations
}
\description{
  Computes the areas of successive morphological dilations.
}
\usage{
  dilated.areas(X, r, W=as.owin(X), ..., constrained=TRUE, exact = FALSE)
}
\arguments{
  \item{X}{
    Object to be dilated.
    A point pattern (object of class \code{"ppp"}),
    a line segment pattern (object of class \code{"psp"}),
    or a window (object of class \code{"owin"}).
  }
  \item{r}{
    Numeric vector of radii for the dilations. 
  }
  \item{W}{
    Window (object of class \code{"owin"}) inside which the areas
    will be computed, if \code{constrained=TRUE}.
  }
  \item{\dots}{Ignored.}
  \item{constrained}{
    Logical flag indicating whether areas should be restricted
    to the window \code{W}.
  }
  \item{exact}{
    Logical flag indicating whether areas should be computed
    using analytic geometry (which is slower but more accurate).
    Currently available only when \code{X} is a point pattern.
  }
}
\details{
  This function computes the areas of the dilations of \code{X}
  by each of the radii \code{r[i]}. Areas may also be computed
  inside a specified window \code{W}.

  The morphological dilation of a set \eqn{X} by a distance \eqn{r > 0}
  is the subset 
  consisting of all points \eqn{x}{x} such that the
  distance from \eqn{x} to \eqn{X} is less than 
  or equal to \eqn{r}.

  When \code{X} is a point pattern, the dilation by a distance
  \eqn{r} is the union of
  discs of radius \eqn{r} centred at the points of \code{X}.

  The argument \code{r} should be a vector of nonnegative numbers.

  If \code{exact=TRUE} and if \code{X} is a point pattern,
  then the areas are computed using analytic geometry, which is
  slower but much more accurate. Otherwise the computation is performed
  using \code{\link{distmap}}.

  To compute the dilated object itself, use \code{\link{dilation}}.
}
\seealso{
  \code{\link{owin}},
  \code{\link{as.owin}},
  \code{\link{dilation}},
  \code{\link{eroded.areas}}
}
\examples{
  X <- runifpoint(10)
  a <- dilated.areas(X, c(0.1,0.2), W=square(1), exact=TRUE)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
