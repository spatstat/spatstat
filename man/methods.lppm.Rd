\name{methods.lppm}
\alias{methods.lppm} %DoNotExport
\alias{coef.lppm}
\alias{extractAIC.lppm}
\alias{formula.lppm}
\alias{logLik.lppm}
\alias{nobs.lppm}
\alias{print.lppm}
\alias{summary.lppm}
\alias{terms.lppm}
\alias{update.lppm}
\alias{vcov.lppm}
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

  \method{print}{lppm}(x, ...)

  \method{summary}{lppm}(object, ...)

  \method{terms}{lppm}(x, ...)

  \method{update}{lppm}(object, ...)

  \method{vcov}{lppm}(object, ...)

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
}
\details{
  These are methods for the generic commands
  \code{\link{coef}},
  \code{\link{extractAIC}},
  \code{\link{formula}},
  \code{\link{logLik}},
  \code{\link{nobs}},
  \code{\link{print}}, 
  \code{\link{summary}}, 
  \code{\link{terms}},
  \code{\link{update}} and
  \code{\link{vcov}}
  for the class \code{"lppm"}. 
}
\value{
  See the default methods.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\seealso{
  \code{\link{lppm}},
  \code{\link{plot.lppm}}.
}
\examples{
  example(lpp)
  fit <- lppm(X, ~x)
  print(fit)
  coef(fit)
  formula(fit)
  terms(fit)
  logLik(fit)
  nobs(fit)
  extractAIC(fit)
  update(fit, ~1)
  vcov(fit)
}
\keyword{spatial}
\keyword{models}
