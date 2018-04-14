\name{measureContinuous}
\alias{measureContinuous}
\alias{measureDiscrete}
\title{
  Discrete and Continuous Components of a Measure
}
\description{
  Given a measure \code{A} (object of class \code{"msr"})
  these functions find the discrete and continuous parts
  of \code{A}.
}
\usage{
measureDiscrete(x)
measureContinuous(x)
}
\arguments{
  \item{x}{
    A measure (object of class \code{"msr"}).
  }
}
\details{
  The functions \code{measureDiscrete} and \code{measureContinuous}
  return the discrete and continuous components, respectively,
  of a measure.

  If \code{x} is a measure, then \code{measureDiscrete(x)} is a measure
  consisting only of the discrete (atomic) component of \code{x},
  and \code{measureContinuous(x)} is a measure consisting only of the
  continuous (diffuse) component of \code{x}.
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
  \code{\link{msr}}, \code{\link{with.msr}}, \code{\link{split.msr}},
  \code{\link{measurePositive}}
}
\examples{
   X <- rpoispp(function(x,y) { exp(3+3*x) })
   fit <- ppm(X, ~x+y)
   rp <- residuals(fit, type="pearson")

   rp
   measureDiscrete(rp)
   measureContinuous(rp)
}
\keyword{spatial}
\keyword{math}
