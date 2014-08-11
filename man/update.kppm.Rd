\name{update.kppm}
\alias{update.kppm}
\title{Update a Fitted Cluster Point Process Model}
\description{
  \code{update} method for class \code{"kppm"}.
}
\usage{
 \method{update}{kppm}(object, trend = ~1, ..., clusters = NULL)
}
\arguments{
  \item{object}{
    Fitted cluster point process model.
    An object of class \code{"kppm"},
    obtained from \code{\link{kppm}}.
  }
  \item{trend}{
    A formula without a left hand side,
    determining the form of the intensity of the model.
  }
  \item{\dots}{
    Other arguments passed to \code{\link{kppm}}.
  }
  \item{clusters}{
    The type of cluster mechanism. A character string.
    See \code{\link{kppm}}.
  }
}
\details{
  \code{object} should be a fitted cluster point process model,
  obtained from the model-fitting function \code{\link{kppm}}.
  The model will be updated according to the new arguments provided.

  The argument \code{trend} determines the formula for the intensity
  in the model. It should be an \R formula without a left hand side.
  It may include the symbols \code{.} and \code{+} or \code{-} to
  specify addition or deletion of terms in the current model formula,
  as shown in the Examples below.

  The model is refitted using \code{\link{kppm}}.
}
\value{
  Another fitted cluster point process model (object of
  class \code{"kppm"}.
}
\seealso{
  \code{\link{kppm}},  \code{\link{plot.kppm}},
  \code{\link{predict.kppm}}, \code{\link{simulate.kppm}},
  \code{\link[spatstat:methods.kppm]{methods.kppm}},
  \code{\link{vcov.kppm}}
}
\examples{
 data(redwood)
 fit <- kppm(redwood, ~1, "Thomas")
 fitx <- update(fit, ~ . + x)
 fitM <- update(fit, clusters="MatClust")
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
