\name{rounding} 
\alias{rounding}
\alias{rounding.default}
\alias{rounding.ppp}
\alias{rounding.pp3}
\alias{rounding.ppx}
\title{
  Detect Numerical Rounding 
}
\description{
  Given a numeric vector, or an object containing numeric spatial
  coordinates, determine whether the values have been rounded
  to a certain number of decimal places.
}
\usage{
rounding(x)

\method{rounding}{default}(x)

\method{rounding}{ppp}(x)

\method{rounding}{pp3}(x)

\method{rounding}{ppx}(x)
}
\arguments{
  \item{x}{
    A numeric vector, or an object containing numeric spatial coordinates.
  }
}
\details{
  For a numeric vector \code{x},
  this function determines whether the values have been rounded
  to a certain number of decimal places.
  \itemize{
    \item 
    If the entries of \code{x} are not all integers, then
    \code{rounding(x)} returns the smallest number of digits \code{d}
    after the decimal point
    such that \code{\link[base]{round}(x, digits=d)} is identical to
    \code{x}.
    For example if \code{rounding(x) = 2} then the entries of
    \code{x} are rounded to 2 decimal places, and are multiples of 0.01.
    \item 
    If all the entries of \code{x} are integers, then
    \code{rounding(x)} returns \code{-d}, where
    \code{d} is the smallest number of digits \emph{before} the decimal point
    such that \code{\link[base]{round}(x, digits=-d)} is identical to
    \code{x}.
    For example if \code{rounding(x) = -3} then the entries of
    \code{x} are multiples of 1000.
    If \code{rounding(x) = 0} then the entries of \code{x} are integers
    but not multiples of 10.
    \item
    If all entries of \code{x} are equal to 0, the rounding is
    not determined, and a value of \code{NULL} is returned.
  }
  
  For a point pattern (object of class \code{"ppp"})
  or similar object \code{x} containing numeric spatial
  coordinates, this procedure is applied to the spatial coordinates.
}
\value{
  An integer.
}
\author{
  \adrian
  
  
  and \rolf
  
}
\seealso{
  \code{\link{round.ppp}}
}
\examples{
   rounding(c(0.1, 0.3, 1.2))
   rounding(c(1940, 1880, 2010))
   rounding(0)
   rounding(cells)
}
\keyword{spatial}
\keyword{math}
