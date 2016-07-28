\name{linearpcfinhom}
\alias{linearpcfinhom}
\title{
  Inhomogeneous Linear Pair Correlation Function
}
\description{
  Computes an estimate of the inhomogeneous linear pair correlation function
  for a point pattern on a linear network.
}
\usage{
linearpcfinhom(X, lambda=NULL, r=NULL, ..., correction="Ang",
               normalise=TRUE, normpower=1)
}
\arguments{
  \item{X}{
    Point pattern on linear network (object of class \code{"lpp"}).
  }
  \item{lambda}{
    Intensity values for the point pattern. Either a numeric vector,
    a \code{function}, a pixel image (object of class \code{"im"}) or
    a fitted point process model (object of class \code{"ppm"}
    or \code{"lppm"}).
  }
  \item{r}{
    Optional. Numeric vector of values of the function argument \eqn{r}.
    There is a sensible default.
  }
  \item{\dots}{
    Arguments passed to \code{\link{density.default}}
    to control the smoothing.
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
  \item{normpower}{
    Integer (usually either 1 or 2).
    Normalisation power. See explanation in \code{\link{linearKinhom}}.
  }
}
\details{
  This command computes the inhomogeneous version of the 
  linear pair correlation function from point pattern data on a linear network.

  If \code{lambda = NULL} the result is equivalent to the
  homogeneous pair correlation function \code{\link{linearpcf}}.
  If \code{lambda} is given, then it is expected to provide estimated values
  of the intensity of the point process at each point of \code{X}. 
  The argument \code{lambda} may be a numeric vector (of length equal to
  the number of points in \code{X}), or a \code{function(x,y)} that will be
  evaluated at the points of \code{X} to yield numeric values, 
  or a pixel image (object of class \code{"im"}) or a fitted point 
  process model (object of class \code{"ppm"} or \code{"lppm"}).

  If \code{lambda} is a fitted point process model,
  the default behaviour is to update the model by re-fitting it to
  the data, before computing the fitted intensity.
  This can be disabled by setting \code{update=FALSE}.
  
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
  \adrian
  
  
}
\references{
  Ang, Q.W. (2010) Statistical methodology for spatial point patterns
  on a linear network. MSc thesis, University of Western Australia.

  Ang, Q.W., Baddeley, A. and Nair, G. (2012)
  Geometrically corrected second-order analysis of 
  events on a linear network, with applications to
  ecology and criminology.
  \emph{Scandinavian Journal of Statistics} \bold{39}, 591--617.

  Okabe, A. and Yamada, I. (2001) The K-function method on a network and
  its computational implementation. \emph{Geographical Analysis}
  \bold{33}, 271-290.
}

\seealso{
  \code{\link{linearpcf}},
  \code{\link{linearKinhom}},
  \code{\link{lpp}}
}
\examples{
  data(simplenet)
  X <- rpoislpp(5, simplenet)
  fit <- lppm(X, ~x)
  K <- linearpcfinhom(X, lambda=fit)
  plot(K)
}
\keyword{spatial}
\keyword{nonparametric}


