\name{methods.kppm}
\alias{methods.kppm} %DoNotExport
\alias{coef.kppm}
\alias{formula.kppm}
\alias{print.kppm}
\alias{terms.kppm}
\alias{labels.kppm}
\title{
  Methods for Cluster Point Process Models
}
\description{
  These are methods for the class \code{"kppm"}. 
}
\usage{
\method{coef}{kppm}(object, \dots)
\method{formula}{kppm}(x, \dots)
\method{print}{kppm}(x, ...)
\method{terms}{kppm}(x, \dots)
\method{labels}{kppm}(object, \dots)
}
\arguments{
  \item{x,object}{
    An object of class \code{"kppm"},
    representing a fitted cluster point process model.
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
}
\details{
  These functions are methods for the generic commands
  \code{\link{coef}},
  \code{\link{formula}},
  \code{\link{print}},
  \code{\link{terms}} and
  \code{\link{labels}}
  for the class \code{"kppm"}.

  An object of class \code{"kppm"} represents a fitted
  cluster point process model. It is obtained from \code{\link{kppm}}.
  
  The method \code{coef.kppm} returns the vector of
  \emph{regression coefficients} of the fitted model.
  It does not return the clustering parameters.
}
\value{
  See the help files for the corresponding generic functions.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
}
\seealso{
  \code{\link{kppm}},  \code{\link{plot.kppm}},
  \code{\link{predict.kppm}}, \code{\link{simulate.kppm}},
  \code{\link{update.kppm}}, \code{\link{vcov.kppm}},
  \code{\link{as.ppm.kppm}}.
}
\examples{
  data(redwood)
  fit <- kppm(redwood ~ x, "MatClust")
  coef(fit)
  formula(fit)
  tf <- terms(fit)
  labels(fit)
}
\keyword{spatial}
\keyword{methods}
