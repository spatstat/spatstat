\name{summary.kppm}
\alias{summary.kppm}
\alias{print.summary.kppm}
\title{Summarizing a Fitted Cox or Cluster Point Process Model}
\description{
  \code{summary} method for class \code{"kppm"}.
}
\usage{
  \method{summary}{kppm}(object, \dots, quick=FALSE)

  \method{print}{summary.kppm}(x, \dots)
}
\arguments{
  \item{object}{
    A fitted Cox or cluster point process model (object of
    class \code{"kppm"}).
  }
  \item{quick}{Logical value controlling the scope of the summary.}
  \item{\dots}{Arguments passed to \code{\link{summary.ppm}} or
    \code{\link{print.summary.ppm}} controlling the treatment of
    the trend component of the model.}
  \item{x}{Object of class \code{"summary.kppm"} as returned by
    \code{summary.kppm}.
  }
}
\details{
  This is a method for the generic \code{\link{summary}}
  for the class \code{"kppm"}. An object of class \code{"kppm"}
  describes a fitted Cox or cluster point process model.
  See \code{\link{kppm}}.
  
  \code{summary.kppm} extracts information about the
  type of model that has been fitted, the data to which the model was
  fitted, and the values of the fitted coefficients.

  \code{print.summary.kppm} prints this information in a
  comprehensible format.

  In normal usage, \code{print.summary.kppm} is invoked implicitly
  when the user calls \code{summary.kppm} without assigning its value
  to anything. See the examples.

  You can also type \code{coef(summary(object))} to extract a table
  of the fitted coefficients of the point process model \code{object}
  together with standard errors and confidence limits. 
}
\value{
  \code{summary.kppm} returns an object of class \code{"summary.kppm"},
  while \code{print.summary.kppm} returns \code{NULL}.
}
\examples{
 fit <- kppm(redwood ~ 1, "Thomas")
 summary(fit)
 coef(summary(fit))
}
\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{methods}
\keyword{models}
