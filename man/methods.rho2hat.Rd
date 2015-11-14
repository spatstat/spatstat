\name{methods.rho2hat}
\alias{methods.rho2hat} %DoNotExport
\alias{predict.rho2hat}
\alias{print.rho2hat}
\alias{plot.rho2hat}
\title{
  Methods for Intensity Functions of Two Spatial Covariates
}
\description{
  These are methods for the class \code{"rho2hat"}. 
}
\usage{
\method{plot}{rho2hat}(x, \dots, do.points=FALSE)

\method{print}{rho2hat}(x, \dots)

\method{predict}{rho2hat}(object, \dots, relative=FALSE)
}
\arguments{
  \item{x,object}{
    An object of class \code{"rho2hat"}.
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
  \item{do.points}{
    Logical value indicating whether to plot the observed values of
    the covariates at the data points.
  }
  \item{relative}{
    Logical value indicating whether to compute the
    estimated point process intensity (\code{relative=FALSE})
    or the relative risk (\code{relative=TRUE}) in the case
    of a relative risk estimate.
  }
}
\details{
  These functions are methods for the generic commands
  \code{\link{print}},
  \code{\link{predict}} and
  \code{\link{plot}} 
  for the class \code{"rho2hat"}.

  An object of class \code{"rho2hat"} is an estimate
  of the intensity of a point process, as a function of two
  given spatial covariates. See \code{\link{rho2hat}}.
  
  The method \code{plot.rho2hat} displays the estimated function
  \eqn{\rho}{rho} using \code{\link{plot.fv}}, and optionally
  adds a \code{\link{rug}} plot of the observed values of the covariate.
  In this plot the two axes represent possible values of the two covariates.

  The method \code{predict.rho2hat} computes a pixel image of the
  intensity \eqn{\rho(Z_1(u), Z_2(u))}{rho(Z1(u), Z2(u))}
  at each spatial location \eqn{u}, where \eqn{Z_1(u)}{Z1(u)}
  and \eqn{Z_2(u)}{Z2(u)} are the two spatial covariates.
}
\value{
  For \code{predict.rho2hat} the value is a pixel image
  (object of class \code{"im"}). 
  For other functions, the value is \code{NULL}.
}
\author{
  \adrian
  
  
}
\seealso{
  \code{\link{rho2hat}}
}
\examples{
  r2 <- with(bei.extra, rho2hat(bei, elev, grad))
  r2
  plot(r2)
  plot(predict(r2))
}
\keyword{spatial}
\keyword{methods}
