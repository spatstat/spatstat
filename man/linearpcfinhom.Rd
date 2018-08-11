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
               normalise=TRUE, normpower=1,
	       update = TRUE, leaveoneout = TRUE,
	       ratio = FALSE)
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
    points, raised to \code{normpower}), which reduces the sampling variability.
    If \code{FALSE}, the denominator is the length of the network.
  }
  \item{normpower}{
    Integer (usually either 1 or 2).
    Normalisation power. See explanation in \code{\link{linearKinhom}}.
  }
  \item{update}{
    Logical value indicating what to do when \code{lambda} is a fitted model
    (class \code{"lppm"} or \code{"ppm"}).
    If \code{update=TRUE} (the default),
    the model will first be refitted to the data \code{X}
    (using \code{\link{update.lppm}} or \code{\link{update.ppm}})
    before the fitted intensity is computed.
    If \code{update=FALSE}, the fitted intensity of the
    model will be computed without re-fitting it to \code{X}.
  }
  \item{leaveoneout}{
    Logical value (passed to \code{\link{fitted.lppm}} or
    \code{\link{fitted.ppm}}) specifying whether to use a
    leave-one-out rule when calculating the intensity,
    when \code{lambda} is a fitted model.
    Supported only when \code{update=TRUE}.
  }
  \item{ratio}{
    Logical. 
    If \code{TRUE}, the numerator and denominator of
    each estimate will also be saved,
    for use in analysing replicated point patterns.
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

  The bandwidth for smoothing the pairwise distances
  is determined by arguments \code{\dots}
  passed to \code{\link{density.default}}, mainly the arguments
  \code{bw} and \code{adjust}. The default is
  to choose the bandwidth by Silverman's rule of thumb 
  \code{bw="nrd0"} explained in \code{\link{density.default}}.
}
\value{
  Function value table (object of class \code{"fv"}).

  If \code{ratio=TRUE} then the return value also has two
  attributes called \code{"numerator"} and \code{"denominator"}
  which are \code{"fv"} objects
  containing the numerators and denominators of each
  estimate of \eqn{g(r)}. 
}
\author{
  Ang Qi Wei \email{aqw07398@hotmail.com} and
  \adrian.
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
  fit <- lppm(X ~x)
  K <- linearpcfinhom(X, lambda=fit)
  plot(K)
}
\keyword{spatial}
\keyword{nonparametric}


