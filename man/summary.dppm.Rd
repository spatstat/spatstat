\name{summary.dppm}
\alias{summary.dppm}
\alias{print.summary.dppm}
\title{Summarizing a Fitted Determinantal Point Process Model}
\description{
  \code{summary} method for class \code{"dppm"}.
}
\usage{
  \method{summary}{dppm}(object, \dots, quick=FALSE)

  \method{print}{summary.dppm}(x, \dots)
}
\arguments{
  \item{object}{
    A fitted determinantal point process model (object of
    class \code{"dppm"}).
  }
  \item{quick}{Logical value controlling the scope of the summary.}
  \item{\dots}{Arguments passed to \code{\link{summary.ppm}} or
    \code{\link{print.summary.ppm}} controlling the treatment of
    the trend component of the model.}
  \item{x}{Object of class \code{"summary.dppm"} as returned by
    \code{summary.dppm}.
  }
}
\details{
  This is a method for the generic \code{\link{summary}}
  for the class \code{"dppm"}. An object of class \code{"dppm"}
  describes a fitted determinantal point process model.
  See \code{\link{dppm}}.
  
  \code{summary.dppm} extracts information about the
  type of model that has been fitted, the data to which the model was
  fitted, and the values of the fitted coefficients.

  \code{print.summary.dppm} prints this information in a
  comprehensible format.

  In normal usage, \code{print.summary.dppm} is invoked implicitly
  when the user calls \code{summary.dppm} without assigning its value
  to anything. See the examples.
}
\value{
  \code{summary.dppm} returns an object of class \code{"summary.dppm"},
  while \code{print.summary.dppm} returns \code{NULL}.

  The result of \code{summary.dppm} includes at least the
  following components:
  \item{Xname}{character string name of the original point pattern data}
  \item{stationary}{logical value indicating whether the model is
    stationary}
  \item{trend}{Object of class \code{summary.ppm} summarising the trend}
  \item{repul}{Repulsiveness index}
}
\examples{
  jpines <- residualspaper$Fig1
  \testonly{
     # smaller dataset for testing
    jpines <- jpines[c(TRUE,FALSE)]
  }
 fit <- dppm(jpines ~ 1, dppGauss)
 summary(fit)
}
\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{methods}
\keyword{models}
