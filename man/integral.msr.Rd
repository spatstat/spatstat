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
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  The integral (total value of the measure over its domain) is calculated.
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
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
