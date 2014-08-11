\name{leverage.ppm}
\alias{leverage}
\alias{leverage.ppm}
\title{
  Leverage Measure for Spatial Point Process Model
}
\description{
  Computes the leverage measure for a fitted spatial point process model.
}
\usage{
leverage(model, ...)
\method{leverage}{ppm}(model, ..., drop = FALSE, iScore=NULL, iHessian=NULL, iArgs=list())
}
\arguments{
  \item{model}{
    Fitted point process model (object of class \code{"ppm"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{drop}{
    Logical. Whether to include (\code{drop=FALSE}) or
    exclude (\code{drop=TRUE}) contributions from quadrature
    points that were not used to fit the model.
  }
  \item{iScore,iHessian}{
    Components of the score vector and Hessian matrix for
    the irregular parameters, if required. See Details.
  }
  \item{iArgs}{
    List of extra arguments for the functions \code{iScore},
    \code{iHessian} if required.
  }
}
\details{
  The function \code{leverage} is generic, and
  \code{leverage.ppm} is the method for objects of class \code{"ppm"}.
  
  Given a fitted spatial point process model \code{model},
  the function \code{leverage.ppm} computes the leverage of the model,
  described in Baddeley, Chang and Song (2011).
  
  The leverage of a spatial point process model
  is a function of spatial location, and is typically
  displayed as a colour pixel image. 
  The leverage value \eqn{h(u)} at a spatial location \eqn{u} represents the
  change in the fitted trend of the fitted point process model that would have
  occurred if a data point were to have occurred at the location \eqn{u}. 
  A relatively large value of \eqn{h()} indicates a 
  part of the space where the data have a \emph{potentially}
  strong effect on the fitted model (specifically, a strong effect
  on the intensity or trend of the fitted model) due to the values
  of the covariates. 
  
  If the point process model trend has irregular parameters that were
  fitted (using \code{\link{ippm}})
  then the leverage calculation requires the first and second
  derivatives of the log trend with respect to the irregular parameters. 
  The argument \code{iScore} should be a list,
  with one entry for each irregular parameter, of \R functions that compute the
  partial derivatives of the log trend (i.e. log intensity or
  log conditional intensity) with respect to each irregular
  parameter. The argument \code{iHessian} should be a list,
  with \eqn{p^2} entries where \eqn{p} is the number of irregular
  parameters, of \R functions that compute the second order
  partial derivatives of the log trend with respect to each
  pair of irregular parameters.  

  The result of \code{leverage.ppm} is an object of
  class \code{"leverage.ppm"}. It can be plotted
  (by \code{\link{plot.leverage.ppm}}) or converted to a pixel
  image by \code{as.im} (see \code{\link{as.im.leverage.ppm}}).
}
\value{
  An object of class \code{"leverage.ppm"} that can be plotted
  by \code{\link{plot.leverage.ppm}}.
}
\references{
  Baddeley, A. and Chang, Y.M. and Song, Y. (2011)
  Leverage and influence diagnostics for spatial point process models.
  \emph{Scandinavian Journal of Statistics}, in press.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\seealso{
  \code{\link{influence.ppm}},
  \code{\link{dfbetas.ppm}},
  \code{\link{plot.leverage.ppm}}
}
\examples{
   X <- rpoispp(function(x,y) { exp(3+3*x) })
   fit <- ppm(X, ~x+y)
   plot(leverage(fit))
}
\keyword{spatial}
\keyword{models}
