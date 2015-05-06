\name{integral.msr}
\alias{integral.msr}
\title{
  Integral of a Measure
}
\description{
  Computes the integral (total value) of a measure over its domain.
}
\usage{
\method{integral}{msr}(f, domain=NULL, \dots)
}
\arguments{
  \item{f}{
    A signed measure or vector-valued measure (object of class
    \code{"msr"}).
  }
  \item{domain}{
    Optional window specifying the domain of integration.
    Alternatively a tessellation.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  The integral (total value of the measure over its domain) is calculated.

  If \code{domain} is a window (class \code{"owin"}) then the integration
  will be restricted to this window. If \code{domain} is a tessellation
  (class \code{"tess"}) then the integral of \code{f} in each
  tile of \code{domain} will be computed.
}
\value{
  A numeric value (for a signed measure)
  or a vector of values (for a vector-valued measure).
}
\seealso{
  \code{\link{msr}},
  \code{\link{integral}}
}
\examples{
   fit <- ppm(cells, ~x)
   rp <- residuals(fit)
   integral(rp)
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{math}
