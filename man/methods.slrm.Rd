\name{methods.slrm}
\alias{methods.slrm}
\alias{formula.slrm}
\alias{print.slrm}
\alias{terms.slrm}
\alias{labels.slrm}
\title{
  Methods for Spatial Logistic Regression Models
}
\description{
  These are methods for the class \code{"slrm"}. 
}
\usage{
\method{formula}{slrm}(x, \dots)
\method{print}{slrm}(x, ...)
\method{terms}{slrm}(x, \dots)
\method{labels}{slrm}(object, \dots)
}
\arguments{
  \item{x,object}{
    An object of class \code{"slrm"},
    representing a fitted cluster point process model.
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
}
\details{
  These functions are methods for the generic commands
  \code{\link{formula}},
  \code{\link{print}},
  \code{\link{terms}} and
  \code{\link{labels}}
  for the class \code{"slrm"}.

  An object of class \code{"slrm"} represents a fitted
  spatial logistic regression model. It is obtained from \code{\link{slrm}}.
}
\value{
  See the help files for the corresponding generic functions.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\seealso{
  \code{\link{slrm}},  \code{\link{plot.slrm}},
  \code{\link{predict.slrm}}, \code{\link{simulate.slrm}},
  \code{\link{update.slrm}}, \code{\link{vcov.slrm}},
  \code{\link{coef.slrm}}.
}
\examples{
  data(redwood)
  fit <- slrm(redwood ~ x)
  coef(fit)
  formula(fit)
  tf <- terms(fit)  
  labels(fit)  
}
\keyword{spatial}
\keyword{methods}
