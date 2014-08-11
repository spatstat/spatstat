\name{linearKinhom}
\alias{linearKinhom}
\title{
  Inhomogeneous Linear K Function
}
\description{
  Computes an estimate of the inhomogeneous linear \eqn{K} function
  for a point pattern on a linear network.
}
\usage{
linearKinhom(X, lambda=NULL, r=NULL, ..., correction="Ang", normalise=TRUE)
}
\arguments{
  \item{X}{
    Point pattern on linear network (object of class \code{"lpp"}).
  }
  \item{lambda}{
    Intensity values for the point pattern. Either a numeric vector,
    a \code{function}, a pixel image
    (object of class \code{"im"} or \code{"linim"}) or
    a fitted point process model (object of class \code{"ppm"}
    or \code{"lppm"}).
  }
  \item{r}{
    Optional. Numeric vector of values of the function argument \eqn{r}.
    There is a sensible default.
  }
  \item{\dots}{
    Ignored.
  }
  \item{correction}{
    Geometry correction.
    Either \code{"none"} or \code{"Ang"}. See Details.
  }
  \item{normalise}{
    Logical. If \code{TRUE} (the default), the denominator of the estimator is 
    data-dependent (equal to the sum of the reciprocal intensities at the data
    points), which reduces the sampling variability.
    If \code{FALSE}, the denominator is the length of the network.
  }
}
\details{
  This command computes the inhomogeneous version of the 
  linear \eqn{K} function from point pattern data on a linear network.

  If \code{lambda = NULL} the result is equivalent to the
  homogeneous \eqn{K} function \code{\link{linearK}}.
  If \code{lambda} is given, then it is expected to provide estimated values
  of the intensity of the point process at each point of \code{X}. 
  The argument \code{lambda} may be a numeric vector (of length equal to
  the number of points in \code{X}), or a \code{function(x,y)} that will be
  evaluated at the points of \code{X} to yield numeric values, 
  or a pixel image (object of class \code{"im"}) or a fitted point 
  process model (object of class \code{"ppm"} or \code{"lppm"}).

  If \code{correction="none"}, the calculations do not include
  any correction for the geometry of the linear network.
  If \code{correction="Ang"}, the pair counts are weighted using
  Ang's correction (Ang, 2010). 

}
\value{
  Function value table (object of class \code{"fv"}).
}
\author{
  Ang Qi Wei \email{aqw07398@hotmail.com} and
  Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\references{
  Ang, Q.W. (2010) Statistical methodology for spatial point patterns
  on a linear network. MSc thesis, University of Western Australia.
  
  Ang, Q.W., Baddeley, A. and Nair, G. (2012)
  Geometrically corrected second-order analysis of 
  events on a linear network, with applications to
  ecology and criminology.
  To appear in \emph{Scandinavian Journal of Statistics}.
}
\seealso{
  \code{\link{lpp}}
}
\examples{
  data(simplenet)
  X <- rpoislpp(5, simplenet)
  fit <- lppm(X, ~x)
  K <- linearKinhom(X, lambda=fit)
  plot(K)
}
\keyword{spatial}
\keyword{nonparametric}


