\name{Kmodel.kppm}
\alias{Kmodel.kppm}
\alias{pcfmodel.kppm}
\title{K Function or Pair Correlation Function of Cluster Model or Cox model}
\description{
  Returns the theoretical \eqn{K} function or the pair correlation function
  of a cluster point process model or Cox point process model.
}
\usage{
   \method{Kmodel}{kppm}(model, \dots)

   \method{pcfmodel}{kppm}(model, \dots)
}
\arguments{
  \item{model}{
    A fitted cluster point process model (object of
    class \code{"kppm"}) typically obtained from
    the model-fitting algorithm \code{\link{kppm}}. 
  }
  \item{\dots}{
    Ignored.
  }
}
\value{
  A \code{function} in the \R language,
  which takes one argument \code{r}.
}
\details{
  For certain types of point process models, it is possible to
  write down a mathematical expression for the \eqn{K} function
  or the pair correlation function of the model. In particular this
  is possible for a fitted cluster point process model 
  (object of class \code{"kppm"} obtained from \code{\link{kppm}}).
  
  The functions \code{\link{Kmodel}} and \code{\link{pcfmodel}} are generic.
  The functions documented here are the methods for the class \code{"kppm"}.
  
  The return value is a \code{function} in the \R language,
  which takes one argument \code{r}.
  Evaluation of this function, on a numeric vector \code{r},
  yields values of the desired \eqn{K} function or pair correlation
  function at these distance values.
}
\seealso{
  \code{\link{Kest}} or \code{\link{pcf}}
  to estimate the \eqn{K} function or pair correlation function
  nonparametrically from data.

  \code{\link{kppm}} to fit cluster models.

  \code{\link{Kmodel}} for the generic functions.

  \code{\link{Kmodel.ppm}} for the method for Gibbs processes.
}
\examples{
  data(redwood)
  fit <- kppm(redwood, ~x, "MatClust")
  K <- Kmodel(fit)
  K(c(0.1, 0.2))
  curve(K(x), from=0, to=0.25)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{models}
