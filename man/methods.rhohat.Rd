\name{methods.rhohat}
\alias{methods.rhohat} %DoNotExport
\alias{print.rhohat}
\alias{plot.rhohat}
\alias{predict.rhohat}
\alias{simulate.rhohat}
\title{
  Methods for Intensity Functions of Spatial Covariate
}
\description{
  These are methods for the class \code{"rhohat"}. 
}
\usage{
\method{print}{rhohat}(x, ...)

\method{plot}{rhohat}(x, ..., do.rug=TRUE)

\method{predict}{rhohat}(object, ..., relative=FALSE)

\method{simulate}{rhohat}(object, nsim=1, ..., drop=TRUE)
}
\arguments{
  \item{x,object}{
    An object of class \code{"rhohat"} representing 
    a smoothed estimate of the intensity function of a point process.
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
  \item{do.rug}{
    Logical value indicating whether to plot the observed values of
    the covariate as a rug plot along the horizontal axis.
  }
  \item{relative}{
    Logical value indicating whether to compute the
    estimated point process intensity (\code{relative=FALSE})
    or the relative risk (\code{relative=TRUE}) in the case
    of a relative risk estimate.
  }
  \item{nsim}{
    Number of simulations to be generated.
  }
  \item{drop}{
    Logical value indicating what to do when \code{nsim=1}.
    If \code{drop=TRUE} (the default), a point pattern is returned.
    If \code{drop=FALSE}, a list of length 1 containing a point pattern
    is returned.
  }
}
\details{
  These functions are methods for the generic commands
  \code{\link{print}},
  \code{\link[graphics]{plot}},
  \code{\link[stats]{predict}} and
  \code{\link[stats]{simulate}}
  for the class \code{"rhohat"}.

  An object of class \code{"rhohat"} is an estimate
  of the intensity of a point process, as a function of a
  given spatial covariate. See \code{\link{rhohat}}.
  
  The method \code{plot.rhohat} displays the estimated function
  \eqn{\rho}{rho} using \code{\link{plot.fv}}, and optionally
  adds a \code{\link{rug}} plot of the observed values of the covariate.
  
  The method \code{predict.rhohat} computes a pixel image of the
  intensity \eqn{\rho(Z(u))}{rho(Z(u))} at each spatial location
  \eqn{u}, where \eqn{Z} is the spatial covariate.

  The method \code{simulate.rhohat} invokes \code{predict.rhohat}
  to determine the predicted intensity, and then simulates a
  Poisson point process with this intensity.
}
\value{
  For \code{predict.rhohat} the value is a pixel image
  (object of class \code{"im"} or \code{"linim"}).
  For \code{simulate.rhohat} the value is a point pattern
  (object of class \code{"ppp"} or \code{"lpp"}).
  For other functions, the value is \code{NULL}.
}
\author{
  \adrian
}
\seealso{
  \code{\link{rhohat}}
}
\examples{
  X <-  rpoispp(function(x,y){exp(3+3*x)})
  rho <- rhohat(X, function(x,y){x})
  rho
  plot(rho)
  Y <- predict(rho)
  plot(Y)
  plot(simulate(rho), add=TRUE)
  # 
  fit <- ppm(X, ~x)
  rho <- rhohat(fit, "y")
  opa <- par(mfrow=c(1,2))
  plot(predict(rho))
  plot(predict(rho, relative=TRUE))
  par(opa)
}
\keyword{spatial}
\keyword{methods}
