\name{formula.ppm}
\alias{formula.ppm}
\alias{terms.ppm}
\title{
  Model Formulae for Gibbs Point Process Models
}
\description{
  Extract the trend formula, or the terms in the trend formula,
  in a fitted Gibbs point process model.
}
\usage{
\method{formula}{ppm}(x, \dots)
\method{terms}{ppm}(x, \dots)
}
\arguments{
  \item{x}{
    An object of class \code{"ppm"},
    representing a fitted point process model.
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
}
\details{
  These functions are methods for the generic commands
  \code{\link{formula}} and 
  \code{\link{terms}} 
  for the class \code{"ppm"}.

  An object of class \code{"ppm"} represents a fitted
  Poisson or Gibbs point process model.
  It is obtained from the model-fitting function \code{\link{ppm}}.

  The method \code{formula.ppm} extracts the trend formula from the
  fitted model \code{x} (the formula originally specified as the
  argument \code{trend} to \code{\link{ppm}}).
  The method \code{terms.ppm} extracts the individual
  terms in the trend formula.
}
\value{
  See the help files for the corresponding generic functions.
}
\author{
  \adrian
  
  
}
\seealso{
  \code{\link{ppm}},  
  \code{\link{as.owin}},
  \code{\link{coef.ppm}},
  \code{\link{extractAIC.ppm}},
  \code{\link{fitted.ppm}},
  \code{\link{logLik.ppm}},
  \code{\link{model.frame.ppm}},
  \code{\link{model.matrix.ppm}},
  \code{\link{plot.ppm}},
  \code{\link{predict.ppm}},
  \code{\link{residuals.ppm}},
  \code{\link{simulate.ppm}},
  \code{\link{summary.ppm}},
  \code{\link{update.ppm}},
  \code{\link{vcov.ppm}}.
}
\examples{
  data(cells)
  fit <- ppm(cells, ~x)
  formula(fit)
  terms(fit)  
}
\keyword{spatial}
\keyword{methods}
