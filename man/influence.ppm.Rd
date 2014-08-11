\name{influence.ppm}
\alias{influence.ppm}
\title{
  Influence Measure for Spatial Point Process Model
}
\description{
  Computes the influence measure for a fitted spatial point process model.
}
\usage{
\method{influence}{ppm}(model, ..., drop = FALSE, iScore=NULL, iHessian=NULL, iArgs=list())
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
  Given a fitted spatial point process model \code{model},
  this function computes the influence measure
  described in Baddeley, Chang and Song (2013).
   
  The function \code{\link[stats]{influence}} is generic,
  and \code{influence.ppm} is the method for objects of class
  \code{"ppm"} representing point process models.

  The influence of a point process model is a value attached to each data point
  (i.e. each point of the point pattern to which the \code{model}
  was fitted).
  The influence value \eqn{s(x_i)}{s(x[i])} at a data point
  \eqn{x_i}{x[i]} represents the change in the maximised log (pseudo)likelihood
  that occurs when the point \eqn{x_i}{x[i]} is deleted.
  A relatively large value of \eqn{s(x_i)}{s(x[i])} indicates a 
  data point with a large influence on the fitted model.
  
  If the point process model trend has irregular parameters that were
  fitted (using \code{\link{ippm}})
  then the influence calculation requires the first and second
  derivatives of the log trend with respect to the irregular parameters. 
  The argument \code{iScore} should be a list,
  with one entry for each irregular parameter, of \R functions that compute the
  partial derivatives of the log trend (i.e. log intensity or
  log conditional intensity) with respect to each irregular
  parameter. The argument \code{iHessian} should be a list,
  with \eqn{p^2} entries where \eqn{p} is the number of irregular
  parameters, of \R functions that compute the second order
  partial derivatives of the
  log trend with respect to each pair of irregular parameters.
  
  The result of \code{influence.ppm} is 
  an object of class \code{"influence.ppm"}. It can be plotted
  (by \code{\link{plot.influence.ppm}}), or converted to a marked
  point pattern by \code{as.ppp} (see \code{\link{as.ppp.influence.ppm}}).
}
\value{
  An object of class \code{"influence.ppm"} that can be plotted
  by \code{\link{plot.influence.ppm}}.
}
\references{
  Baddeley, A. and Chang, Y.M. and Song, Y. (2013)
  Leverage and influence diagnostics for spatial point process models.
  \emph{Scandinavian Journal of Statistics} \bold{40}, 86--104.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\seealso{
  \code{\link{leverage.ppm}},
  \code{\link{dfbetas.ppm}},
  \code{\link{plot.influence.ppm}}
}
\examples{
   X <- rpoispp(function(x,y) { exp(3+3*x) })
   fit <- ppm(X, ~x+y)
   plot(influence(fit))
}
\keyword{spatial}
\keyword{models}
