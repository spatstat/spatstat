\name{measureVariation}
\alias{measureVariation}
\alias{measurePositive}
\alias{measureNegative}
\title{
  Positive and Negative Parts, and Variation, of a Measure
}
\description{
  Given a measure \code{A} (object of class \code{"msr"})
  these functions find the positive part, negative part and variation
  of \code{A}.
}
\usage{
measurePositive(x)
measureNegative(x)
measureVariation(x)
}
\arguments{
  \item{x}{
    A measure (object of class \code{"msr"}).
  }
}
\details{
  The functions \code{measurePositive} and \code{measureNegative}
  return the positive and negative parts of the measure,
  and \code{measureVariation} returns the variation.

  If \eqn{\mu} is a signed measure,
  it can be represented as
  \deqn{\mu = \mu_{+} - \mu_{-}}{\mu = \mu[+] - \mu[-]}
  where \eqn{\mu_{+}}{\mu[+]} and \eqn{\mu_{-}}{\mu[-]}
  are \emph{nonnegative} measures called the positive and negative
  parts.
  In a nutshell, the positive part of \eqn{\mu}
  consists of all positive contributions or increments,
  and the negative part consists of all negative contributions.

  The variation \eqn{|\mu|} is defined by 
  \deqn{\mu = \mu_{+} + \mu_{-}}{\mu = \mu[+] + \mu[-]}
  and is also a nonnegative measure.
}
\value{
  Another measure (object of class \code{"msr"})
  on the same spatial domain.
}
\references{
  Halmos, P.R. (1950) \emph{Measure Theory}. Van Nostrand.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{msr}}, \code{\link{with.msr}}, \code{\link{split.msr}}
}
\examples{
   X <- rpoispp(function(x,y) { exp(3+3*x) })
   fit <- ppm(X, ~x+y)
   rp <- residuals(fit, type="pearson")

   measurePositive(rp)
   measureNegative(rp)
   measureVariation(rp)

   # total variation norm
   integral(measureVariation(rp))
}
\keyword{spatial}
\keyword{math}
