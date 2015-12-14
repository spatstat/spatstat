\name{rho2hat}
\alias{rho2hat}
\title{
  Smoothed Relative Density of Pairs of Covariate Values
}
\description{
  Given a point pattern and two spatial covariates \eqn{Z_1}{Z1} and
  \eqn{Z_2}{Z2}, construct a smooth estimate of the relative risk of
  the pair \eqn{(Z_1,Z_2)}{(Z1, Z2)}.
}
\usage{
rho2hat(object, cov1, cov2, ..., method=c("ratio", "reweight"))
}
\arguments{
  \item{object}{
    A point pattern (object of class \code{"ppp"}),
    a quadrature scheme (object of class \code{"quad"})
    or a fitted point process model (object of class \code{"ppm"}).
  }
  \item{cov1,cov2}{
    The two covariates.
    Each argument is either a \code{function(x,y)} or a pixel image (object of
    class \code{"im"}) providing the values of the covariate at any
    location, or one of the strings \code{"x"} or \code{"y"}
    signifying the Cartesian coordinates.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{density.ppp}} to smooth
    the scatterplots.
  }
  \item{method}{
    Character string determining the smoothing method. See Details.
  }
}
\details{
  This is a bivariate version of \code{\link{rhohat}}.
  
  If \code{object} is a point pattern, this command
  produces a smoothed version of the scatterplot of
  the values of the covariates \code{cov1} and \code{cov2}
  observed at the points of the point pattern. 

  The covariates \code{cov1,cov2} must have continuous values.
  
  If \code{object} is a fitted point process model, suppose \code{X} is
  the original data point pattern to which the model was fitted. Then
  this command assumes \code{X} is a realisation of a Poisson point
  process with intensity function of the form
  \deqn{
    \lambda(u) = \rho(Z_1(u), Z_2(u)) \kappa(u)
  }{
    lambda(u) = rho(Z1(u), Z2(u)) * kappa(u)
  }
  where \eqn{\kappa(u)}{kappa(u)} is the intensity of the fitted model
  \code{object}, and \eqn{\rho(z_1,z_2)}{rho(z1, z2)} is a function
  to be estimated. The algorithm computes a smooth estimate of the
  function \eqn{\rho}{rho}.

  The \code{method} determines how the density estimates will be
  combined to obtain an estimate of \eqn{\rho(z_1, z_2)}{rho(z1, z2)}:
  \itemize{
    \item
    If \code{method="ratio"}, then \eqn{\rho(z_1, z_2)}{rho(z1,z2)} is
    estimated by the ratio of two density estimates.
    The numerator is a (rescaled) density estimate obtained by
    smoothing the points \eqn{(Z_1(y_i), Z_2(y_i))}{(Z1(y[i]), Z2(y[i]))} 
    obtained by evaluating the two covariate \eqn{Z_1, Z_2}{Z1, Z2}
    at the data points \eqn{y_i}{y[i]}. The denominator
    is a density estimate of the reference distribution of
    \eqn{(Z_1,Z_2)}{(Z1, Z2)}.
    \item
    If \code{method="reweight"}, then \eqn{\rho(z_1, z_2)}{rho(z1,z2)} is
    estimated by applying density estimation to the 
    points \eqn{(Z_1(y_i), Z_2(y_i))}{(Z1(y[i]), Z2(y[i]))}
    obtained by evaluating the two covariate \eqn{Z_1, Z_2}{Z1, Z2}
    at the data points \eqn{y_i}{y[i]},
    with weights inversely proportional to the reference density of
    \eqn{(Z_1,Z_2)}{(Z1, Z2)}.
    }
}
\value{
  A pixel image (object of class \code{"im"}). Also
  belongs to the special class \code{"rho2hat"} which has a plot method.
}
\references{
  Baddeley, A., Chang, Y.-M., Song, Y. and Turner, R. (2012)
  Nonparametric estimation of the dependence of a point
  process on spatial covariates.
  \emph{Statistics and Its Interface} \bold{5} (2), 221--236.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
}
\seealso{
  \code{\link{rhohat}},
  \code{\link{methods.rho2hat}}
}
\examples{
  data(bei)
  attach(bei.extra)
  plot(rho2hat(bei, elev, grad))
  fit <- ppm(bei, ~elev, covariates=bei.extra)
  \dontrun{
  plot(rho2hat(fit, elev, grad))
  }
  plot(rho2hat(fit, elev, grad, method="reweight"))
}
\keyword{spatial}
\keyword{models}
