\name{logLik.mppm}
\alias{logLik.mppm}
\title{Log Likelihood for Poisson Point Process Model}
\description{
  Extracts the log likelihood of a Poisson point process model
  that has been fitted to multiple point patterns
}
\usage{
 \method{logLik}{mppm}(object, ...)
}
\arguments{
  \item{object}{Fitted point process model (fitted to multiple point
    patterns). An object of class \code{"mppm"}.
  }
  \item{\dots}{Ignored.}
}
\details{
  The maximised value of the log likelihood for the fitted model
  (as approximated by quadrature using the Berman-Turner approximation)
  is extracted.

  If \code{object} is not a Poisson process, the maximised log
  \emph{pseudolikelihood} is returned, with a warning.
}
\value{
  A numerical value.
}
\seealso{
  \code{\link{mppm}}
}
\author{Adrian Baddeley
  \email{adrian.baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\examples{
  data(waterstriders)
  fit <- mppm(Bugs ~ x, hyperframe(Bugs=waterstriders))
  logLik(fit)
  AIC(fit)
}
\keyword{spatial}
\keyword{models}

