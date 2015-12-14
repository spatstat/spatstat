\name{spatialcdf}
\alias{spatialcdf}
\title{
  Spatial Cumulative Distribution Function
}
\description{
  Compute the spatial cumulative distribution function of a
  spatial covariate, optionally using spatially-varying weights.
}
\usage{
spatialcdf(Z, weights = NULL, normalise = FALSE, ..., W = NULL, Zname = NULL)
}
\arguments{
  \item{Z}{
    Spatial covariate.
    A pixel image or a \code{function(x,y,...)}
  }
  \item{weights}{
    Spatial weighting for different locations.
    A pixel image, a \code{function(x,y,...)}, a window, a constant value,
    or a fitted point process model (object of class \code{"ppm"} or
    \code{"kppm"}).
  }
  \item{normalise}{
    Logical. Whether the weights should be normalised so that they
    sum to 1.
  }
  \item{\dots}{
    Arguments passed to \code{\link{as.mask}} to determine the pixel
    resolution, or extra arguments passed to \code{Z} if it is a function.
  }
  \item{W}{
    Optional window (object of class \code{"owin"}) defining the spatial
    domain.
  }
  \item{Zname}{
    Optional character string for the name of the covariate \code{Z}
    used in plots.
  }
}
\details{
  If \code{weights} is missing or \code{NULL}, it defaults to 1.
  The values of the covariate \code{Z}
  are computed on a grid of pixels. The weighted cumulative distribution
  function of \code{Z} values is computed, taking each value with weight
  equal to the pixel area. The resulting function \eqn{F} is such that
  \eqn{F(t)} is the area of the region of space where
  \eqn{Z \le t}{Z <= t}.

  If \code{weights} is a pixel image or a function, then the
  values of \code{weights} and of the covariate \code{Z}
  are computed on a grid of pixels. The
  \code{weights} are multiplied by the pixel area.
  Then the weighted empirical cumulative distribution function
  of \code{Z} values
  is computed using \code{\link{ewcdf}}. The resulting function
  \eqn{F} is such that \eqn{F(t)} is the total weight (or weighted area)
  of the region of space where \eqn{Z \le t}{Z <= t}.

  If \code{weights} is a fitted point process model, then it should
  be a Poisson process. The fitted intensity of the model,
  and the value of the covariate \code{Z}, are evaluated at the
  quadrature points used to fit the model. The \code{weights} are
  multiplied by the weights of the quadrature points.
  Then the weighted empirical cumulative distribution of \code{Z} values
  is computed using \code{\link{ewcdf}}. The resulting function
  \eqn{F} is such that \eqn{F(t)} is the expected number of points
  in the point process that will fall in the region of space
  where \eqn{Z \le t}{Z <= t}.
  
  If \code{normalise=TRUE}, the function is normalised so that its
  maximum value equals 1, so that it gives the cumulative
  \emph{fraction} of weight or cumulative fraction of points.

  The result can be printed, plotted, and used as a function.
}
\value{
  A cumulative distribution function object
  belonging to the classes \code{"spatialcdf"},
  \code{"ewcdf"}, \code{"ecdf"} and \code{"stepfun"}.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  ,
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}
  
}
\seealso{
  \code{\link{ewcdf}},
  \code{\link{cdf.test}}
}
\examples{
   with(bei.extra, {
     plot(spatialcdf(grad))
     fit <- ppm(bei ~ elev)
     plot(spatialcdf(grad, predict(fit)))
     plot(A <- spatialcdf(grad, fit))
     A(0.1)
  })
}
\keyword{spatial}
\keyword{nonparametric}
