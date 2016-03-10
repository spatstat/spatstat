\name{update.kppm}
\alias{update.kppm}
\title{Update a Fitted Cluster Point Process Model}
\description{
  \code{update} method for class \code{"kppm"}.
}
\usage{
 \method{update}{kppm}(object, \dots)
}
\arguments{
  \item{object}{
    Fitted cluster point process model.
    An object of class \code{"kppm"},
    obtained from \code{\link{kppm}}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{kppm}}.
  }
}
\details{
  \code{object} should be a fitted cluster point process model,
  obtained from the model-fitting function \code{\link{kppm}}.
  The model will be updated according to the new arguments provided.

  If the argument \code{trend} is provided, it determines the
  intensity in the updated model. It should be an \R formula
  (with or without a left hand side). It may include the symbols
  \code{+} or \code{-} to specify addition or deletion of terms
  in the current model formula, as shown in the Examples below.
  The symbol \code{.} refers to the current contents of the
  formula.

  The intensity in the updated model is determined by the
  argument \code{trend} if it is provided, or otherwise by any unnamed
  argument that is a formula, or otherwise by the formula of the
  original model, \code{formula(object)}.

  The spatial point pattern data to which the new model is fitted
  is determined by the left hand side of the updated model formula,
  if this is present. Otherwise it is determined by the argument
  \code{X} if it is provided, or otherwise by any unnamed argument
  that is a point pattern or a quadrature scheme.

  The model is refitted using \code{\link{kppm}}.
}
\value{
  Another fitted cluster point process model (object of
  class \code{"kppm"}.
}
\seealso{
  \code{\link{kppm}},  \code{\link{plot.kppm}},
  \code{\link{predict.kppm}}, \code{\link{simulate.kppm}},
  \code{\link{methods.kppm}},
  \code{\link{vcov.kppm}}
}
\examples{
 fit <- kppm(redwood ~1, "Thomas")
 fitx <- update(fit, ~ . + x)
 fitM <- update(fit, clusters="MatClust")
 fitC <- update(fit, cells)
 fitCx <- update(fit, cells ~ x)
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}
  
}
\keyword{spatial}
\keyword{models}
