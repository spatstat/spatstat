\name{rhohatMonotone}
\alias{rhohatMonotone}
\concept{Resource Selection Function}
\concept{Prospectivity}
\title{
  Monotone Estimate of Intensity as Function of a Covariate
}
\description{
  Computes an estimate of the intensity of a point process,
  as a function of a (continuous) spatial covariate,
  assuming the function is monotone increasing or
  monotone decreasing.
}
\usage{
rhohatMonotone(X, Z, increasing = FALSE)
}
\arguments{
  \item{X}{
    Point pattern dataset (object of class \code{"ppp"}
    or \code{"lpp"}).
  }
  \item{Z}{
    Covariate.
    Either a \code{function(x,y)} or a pixel image (object of
    class \code{"im"}) providing the values of the covariate at any
    location.
    Alternatively one of the strings \code{"x"} or \code{"y"}
    signifying the Cartesian coordinates.
  }
  \item{increasing}{
    Logical value indicating whether the intensity is
    assumed to be increasing as a function of \code{Z}
    (\code{increasing=TRUE}) or decreasing as a function of \code{Z}
    (\code{increasing=FALSE}, the default).
  }
}
\details{
  This command estimates the relationship between
  point process intensity and a given spatial covariate.
  Such a relationship is sometimes called a
  \emph{resource selection function} (if the points are organisms
  and the covariate is a descriptor of habitat) or
  a \emph{prospectivity index} (if the points are mineral deposits
  and the covariate is a geological variable).

  The intensity \eqn{\lambda(u)}{lambda(u)}
  of points at a spatial location \eqn{u}
  is assumed to satisfy
  \deqn{
    \lambda(u) = \rho(Z(u))
  }{
    lambda(u) = rho(Z(u))
  }
  where \eqn{\rho(z)}{rho(z)} is the resource selection function
  or prospectivity index. The function \eqn{\rho(z)}{rho(z)} will be
  estimated.

  This particular command assumes that the intensity
  is a \emph{decreasing} function of \code{Z}
  (or an \emph{increasing} function of \code{Z}, if
  \code{increasing=TRUE}).
  It computes the nonparametric maximum likelihood estimate
  of the function \eqn{\rho(z)}{rho(z)} using the
  maximum upper sets algorithm (Sager, 1982; Baddeley, 2018).

  The command \code{\link{rhohat}} estimates the same function
  \eqn{\rho(z)}{rho(z)}, but using kernel smoothing instead of
  nonparametric maximum likelihood. 
}
\value{
  A step function (object of class \code{"stepfun"}).
}
\references{
  Sager, T.W. (1982) 
  Nonparametric maximum likelihood estimation of
  spatial patterns. \emph{Annals of Statistics} \bold{10}, 1125--1136.

  Baddeley, A. (2018)
  A statistical commentary on mineral prospectivity analysis.
  Chapter 2, pages 25--65
  in \emph{Handbook of Mathematical Geosciences: Fifty Years of IAMG},
  edited by B.S. Daya Sagar, Q. Cheng and F.P. Agterberg.
  Springer, Berlin.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{rhohat}} for a nonparametric smoothing estimator
  of the same function \eqn{\rho}{rho}.
}
\examples{
  rh <- with(murchison, rhohatMonotone(gold, distfun(faults)))
  plot(rh,
       xlab="distance to nearest fault [m]",
       ylab="intensity [deposits per m^2]",
       do.points=FALSE)
}
\keyword{spatial}
\keyword{models}
\keyword{nonparametric}
