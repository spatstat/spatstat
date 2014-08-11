\name{methods.lppm}
\alias{methods.lppm}
\alias{coef.lppm}
\alias{extractAIC.lppm}
\alias{formula.lppm}
\alias{logLik.lppm}
\alias{nobs.lppm}
\alias{print.lppm}
\alias{plot.lppm}
\alias{terms.lppm}
\alias{update.lppm}
\alias{as.linnet.lppm}
\title{
  Methods for Fitted Point Process Models on a Linear Network
}
\description{
  These are methods for the class \code{"lppm"} of fitted point process
  models on a linear network.
}
\usage{
  \method{coef}{lppm}(object, ...)

  \method{extractAIC}{lppm}(fit, ...)

  \method{formula}{lppm}(x, ...)

  \method{logLik}{lppm}(object, ...)

  \method{nobs}{lppm}(object, ...)

  \method{plot}{lppm}(x, ..., type="trend")

  \method{print}{lppm}(x, ...)

  \method{terms}{lppm}(x, ...)

  \method{update}{lppm}(object, ...)

  \method{as.linnet}{lppm}(X, ...)
}
\arguments{
  \item{object,fit,x,X}{
    An object of class \code{"lppm"} representing a fitted point process
    model on a linear network.
  }
  \item{\dots}{
    Arguments passed to other methods, usually the
    method for the class \code{"ppm"}.
  }
  \item{type}{
    Character string (either \code{"trend"} or \code{"cif"})
    determining whether to plot the fitted first order trend
    or the conditional intensity.
  }
}
\details{
  These are methods for the generic commands
  \code{\link{coef}},
  \code{\link{extractAIC}},
  \code{\link{formula}},
  \code{\link{logLik}},
  \code{\link{nobs}},
  \code{\link{plot}},
  \code{\link{print}}, 
  \code{\link{terms}} and
  \code{\link{update}}
  for the class \code{"lppm"}. 
}
\value{
  See the default methods.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\seealso{
  \code{\link{lppm}}
}
\examples{
  example(lpp)
  fit <- lppm(X, ~x)
  print(fit)
  plot(fit)
  coef(fit)
  formula(fit)
  terms(fit)
  logLik(fit)
  nobs(fit)
  extractAIC(fit)
  update(fit, ~1)
}
\keyword{spatial}
\keyword{models}
