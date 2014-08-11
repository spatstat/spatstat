\name{Kmodel}
\alias{Kmodel}
\alias{pcfmodel}
\alias{Kmodel.kppm}
\alias{pcfmodel.kppm}
\title{K function of a model}
\description{
  Returns the theoretical \eqn{K} function or the pair correlation function
  of a point process model.
}
\usage{
   Kmodel(model, \dots)
   pcfmodel(model, \dots)
   \method{Kmodel}{kppm}(model, \dots)
   \method{pcfmodel}{kppm}(model, \dots)
}
\arguments{
  \item{model}{
    A fitted cluster point process model, typically obtained from
    the model-fitting algorithm \code{\link{kppm}}. An object of
    class \code{"kppm"}.
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
  is possible for a fitted cluster point process model (object of
  class \code{"kppm"} obtained from \code{\link{kppm}}).
  
  The functions \code{Kmodel} and \code{pcfmodel} are generic.
  Currently the only method is for the class \code{"kppm"}.
  
  The return value is a \code{function} in the \R language,
  which takes one argument \code{r}.
  Evaluation of this function, on a numeric vector \code{r},
  yields values of the desired \eqn{K} function or pair correlation
  function at these distance values.
}
\seealso{
  \code{\link{kppm}},
  \code{\link{Kest}},
  \code{\link{pcf}}
}
\examples{
  data(redwood)
  fit <- kppm(redwood, ~x, "MatClust")
  K <- Kmodel(fit)
  K(c(0.1, 0.2))
  curve(K(x), from=0, to=0.25)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
