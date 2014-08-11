\name{print.quad}
\alias{print.quad}
\title{Print a Quadrature Scheme}
\description{
  \code{print} method for a quadrature scheme.
}
\usage{
 \method{print}{quad}(x,\dots)
}
\arguments{
  \item{x}{
    A quadrature scheme object, typically obtained from
    \code{\link{quadscheme}}.
    An object of class \code{"quad"}.
  }
  \item{\dots}{Ignored.}
}
\value{
  none.
}
\details{
  This is the \code{print} method for the class \code{"quad"}.
  It prints simple information about the quadrature scheme.

  See \code{\link{quad.object}} for details of the class \code{"quad"}.
}
\seealso{
  \code{\link{quadscheme}},
  \code{\link{quad.object}},
  \code{\link{plot.quad}},
  \code{\link{summary.quad}}
}
\examples{
 data(cells)
 Q <- quadscheme(cells)
 Q
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{print}

