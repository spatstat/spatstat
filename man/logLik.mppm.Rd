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
\references{
  Baddeley, A., Rubak, E. and Turner, R. (2015)
  \emph{Spatial Point Patterns: Methodology and Applications with R}.
  London: Chapman and Hall/CRC Press. 
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}, Ida-Maria Sintorn and Leanne Bischoff.
  Implemented by 
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}
  
}
\examples{
  fit <- mppm(Bugs ~ x, hyperframe(Bugs=waterstriders))
  logLik(fit)
  AIC(fit)
}
\keyword{spatial}
\keyword{models}

