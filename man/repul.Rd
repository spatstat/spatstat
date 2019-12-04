\name{repul.dppm}
\alias{repul}
\alias{repul.dppm}
\title{
  Repulsiveness Index of a Determinantal Point Process Model
}
\description{
  Computes a measure of the degree of repulsion between points
  in a determinantal point process model.
}
\usage{
repul(model, \dots)

\method{repul}{dppm}(model, \dots)
}
\arguments{
  \item{model}{
    A fitted point process model
    of determinantal type (object of class \code{"dppm"}).
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  The repulsiveness index \eqn{\mu}{mu} of a determinantal
  point process model was defined by
  Lavancier, \Moller and Rubak (2015) as
  \deqn{
    \mu = \lambda \int (1- g(x)) \, dx
  }{
    mu = lambda * integral[(1-g(x))]
  }
  where \eqn{\lambda}{\lambda} is the intensity of the model and
  \eqn{g(x)} is the pair correlation function, and
  the integral is taken over all two-dimensional vectors \eqn{x}.

  Values of \eqn{\mu} are dimensionless.
  Larger values of \eqn{\mu}{mu} indicate stronger repulsion
  between points.

  If the model is stationary, the result is a single number.

  If the model is not stationary, 
  the result is a pixel image (obtained by multiplying
  the spatially-varying intensity by the integral defined above).
}
\value{
  A numeric value or a pixel image.
}
\references{
  Lavancier, F., \Moller, J.  and Rubak, E. (2015),
  Determinantal point process models and statistical inference.
  \emph{Journal of Royal Statistical Society:
    Series B (Statistical Methodology)},
  \bold{77}, 853--877.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{dppm}}
}
\examples{
  jpines <- residualspaper$Fig1
  \testonly{
     # smaller dataset for testing
    jpines <- jpines[c(TRUE,FALSE)]
  }
  fit <- dppm(jpines ~ 1, dppGauss)
  repul(fit)
}
\keyword{spatial}
\keyword{models}
