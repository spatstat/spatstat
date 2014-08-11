\name{predict.kppm}
\alias{predict.kppm}
\title{Prediction from a Fitted Cluster Point Process Model}
\description{
  Given a fitted cluster point process model,
  this function computes the fitted intensity.
}
\usage{
  \method{predict}{kppm}(object, ...)
}
\arguments{
  \item{object}{
    Fitted cluster point process model.
    An object of class \code{"kppm"}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{predict.ppm}}.
  }
}
\details{
  This is a method for the generic function \code{\link{predict}}.
  The argument \code{object} should be a cluster point process model
  (object of class \code{"kppm"}) obtained using
  the function \code{\link{kppm}}.
  
  The \emph{intensity} of the fitted model
  is computed, using \code{\link{predict.ppm}}.
}
\value{
  Usually a pixel image
  (object of class \code{"im"}), but see \code{\link{predict.ppm}}
  for details.
}
\seealso{
  \code{\link{kppm}},
  \code{\link{plot.kppm}},
  \code{\link{vcov.kppm}},
  \code{\link{predict.ppm}}
}
\examples{
  data(redwood)
  fit <- kppm(redwood, ~x, "Thomas")
  predict(fit)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
