\name{dfbetas.ppm}
\alias{dfbetas.ppm}
\title{
  Parameter influence measure
}
\description{
  Computes the deletion influence measure for each parameter
  in a fitted point process model.
}
\usage{
\method{dfbetas}{ppm}(model, ..., drop = FALSE, iScore=NULL,
iHessian=NULL, iArgs=NULL)
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
  Given a fitted spatial point process \code{model},
  this function computes the influence measure for each parameter,
  as described in Baddeley, Chang and Song (2013).
  
  This is a method for the generic function \code{\link[stats]{dfbetas}}.

  The influence measure for each parameter \eqn{\theta}{\theta} is a
  signed measure in two-dimensional space. It consists of a discrete
  mass on each data point (i.e. each point in the point pattern to which
  the \code{model} was originally fitted) and a continuous density at
  all locations. The mass at a data point represents the change in the 
  fitted value of the parameter \eqn{\theta}{\theta} that would occur
  if this data point were to be deleted.
  The density at other non-data locations represents the
  effect (on the fitted value of \eqn{\theta}{\theta})
  of deleting these locations (and their associated covariate values)
  from the input to the fitting procedure.
  
  If the point process model trend has irregular parameters that were
  fitted (using \code{\link{ippm}})
  then the influence calculation requires the first and second
  derivatives of the log trend with respect to the irregular parameters. 
  The argument \code{iScore} should be a list,
  with one entry for each irregular parameter,  of \R functions that compute the
  partial derivatives of the log trend (i.e. log intensity or
  log conditional intensity) with respect to each irregular
  parameter. The argument \code{iHessian} should be a list,
  with \eqn{p^2} entries where \eqn{p} is the number of irregular
  parameters, of \R functions that compute the second order
  partial derivatives of the
  log trend with respect to each pair of irregular parameters.  
}
\value{
  An object of class \code{"msr"} representing a signed or vector-valued
  measure. 
}
\references{
  Baddeley, A. and Chang, Y.M. and Song, Y. (2013)
  Leverage and influence diagnostics for spatial point process models.
  \emph{Scandinavian Journal of Statistics} \bold{40}, 86--104.
}
\author{
  \adrian
}
\seealso{
  \code{\link{leverage.ppm}},
  \code{\link{influence.ppm}},
  \code{\link{ppmInfluence}}
}
\examples{
   \testonly{op <- spatstat.options(npixel=32)}
   X <- rpoispp(function(x,y) { exp(3+3*x) })
   fit <- ppm(X, ~x+y)
   \testonly{fit <- ppm(X, ~x+y, nd=16)}
   plot(dfbetas(fit))
   plot(Smooth(dfbetas(fit)))
   \testonly{spatstat.options(op)}
}
\keyword{spatial}
\keyword{models}
