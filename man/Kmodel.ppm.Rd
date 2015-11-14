\name{Kmodel.ppm}
\alias{Kmodel.ppm}
\alias{pcfmodel.ppm}
\title{K Function or Pair Correlation Function of Gibbs Point Process model}
\description{
  Returns the theoretical \eqn{K} function or the pair correlation function
  of a fitted Gibbs point process model.
}
\usage{
   \method{Kmodel}{ppm}(model, \dots)

   \method{pcfmodel}{ppm}(model, \dots)
}
\arguments{
  \item{model}{
    A fitted Poisson or Gibbs point process model (object of
    class \code{"ppm"}) typically obtained from
    the model-fitting algorithm \code{\link{ppm}}. 
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
  This function computes an \emph{approximation} to the \eqn{K} function
  or the pair correlation function of a Gibbs point process.
  
  The functions \code{\link{Kmodel}} and \code{\link{pcfmodel}} are generic.
  The functions documented here are the methods for the class
  \code{"ppm"}.
  
  The approximation is only available for stationary 
  pairwise-interaction models.
  It uses the second order Poisson-saddlepoint approximation
  (Baddeley and Nair, 2012b) which is a combination of
  the Poisson-Boltzmann-Emden and Percus-Yevick approximations.

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

  \code{\link{ppm}} to fit Gibbs models.

  \code{\link{Kmodel}} for the generic functions.

  \code{\link{Kmodel.kppm}} for the method for cluster/Cox processes.
}
\examples{
  fit <- ppm(swedishpines, ~1, Strauss(8))
  p <- pcfmodel(fit)
  K <- Kmodel(fit)
  p(6)
  K(8)
  curve(K(x), from=0, to=15)
}
\references{
  Baddeley, A. and Nair, G. (2012a) 
  Fast approximation of the intensity of Gibbs point processes.
  \emph{Electronic Journal of Statistics} \bold{6} 1155--1169.
  
  Baddeley, A. and Nair, G. (2012b)
  Approximating the moments of a spatial point process.
  \emph{Stat} \bold{1}, 1, 18--30.
  doi: 10.1002/sta4.5
}
\author{\adrian
  
  
  and Gopalan Nair.
}
\keyword{spatial}
\keyword{models}
