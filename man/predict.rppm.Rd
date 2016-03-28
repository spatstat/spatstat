\name{predict.rppm}
\alias{fitted.rppm}
\alias{predict.rppm}
\title{
  Make Predictions From a Recursively Partitioned Point Process Model
}
\description{
  Given a model which has been fitted to point pattern data
  by recursive partitioning, compute the predicted intensity
  of the model.
}
\usage{
\method{predict}{rppm}(object, \dots)

\method{fitted}{rppm}(object, \dots)
}
\arguments{
  \item{object}{
    Fitted point process model of class \code{"rppm"}
    produced by the function \code{\link{rppm}}.
  }
  \item{\dots}{
    Optional arguments passed to \code{\link{predict.ppm}}
    to specify the locations where prediction is required.
    (Ignored by \code{fitted.rppm})
  }
}
\details{
  These functions are methods for the generic functions
  \code{\link[stats]{fitted}} and \code{\link[stats]{predict}}.
  They compute the fitted intensity of a point process model.
  The argument \code{object} should be a fitted point process model
  of class \code{"rppm"} produced by the function \code{\link{rppm}}.

  The \code{fitted} method computes the fitted intensity at the original data
  points, yielding a numeric vector with one entry for each data point.

  The \code{predict} method computes the fitted intensity at
  any locations. By default, predictions are
  calculated at a regular grid of spatial locations, and the result
  is a pixel image giving the predicted intensity values at these
  locations.

  Alternatively, predictions can be performed at other
  locations, or a finer grid of locations, or only at certain specified
  locations, using additional arguments \code{\dots}
  which will be interpreted by \code{\link{predict.ppm}}.
  Common arguments are \code{ngrid} to increase the grid resolution,
  \code{window} to specify the prediction region, and \code{locations}
  to specify the exact locations of predictions.
  See \code{\link{predict.ppm}} for details of these arguments.

  Predictions are computed by evaluating the explanatory covariates at
  each desired location, and applying the recursive partitioning rule
  to each set of covariate values.
}
\value{
  The result of \code{fitted.rppm} is a numeric vector.

  The result of \code{predict.rppm} is a pixel image, a list of pixel images,
  or a numeric vector.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au},
 Rolf Turner \email{r.turner@auckland.ac.nz}
 and Ege Rubak \email{rubak@math.aau.dk}.
}
\seealso{
  \code{\link{rppm}},
  \code{\link{plot.rppm}}  
}
\examples{
    fit <- rppm(unmark(gorillas) ~ vegetation, data=gorillas.extra)
    plot(predict(fit))
    lambdaX <- fitted(fit)
    lambdaX[1:5]
    # Mondriaan pictures
    plot(predict(rppm(redwoodfull ~ x + y)))
    points(redwoodfull)
}
\keyword{spatial}
\keyword{models}
