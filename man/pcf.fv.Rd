\name{pcf.fv}
\alias{pcf.fv}
\title{Pair Correlation Function obtained from K Function}
\description{
  Estimates the pair correlation function of
  a point pattern, given an estimate of the K function.
}
\usage{
  \method{pcf}{fv}(X, \dots, method="c")
}
\arguments{
  \item{X}{
    An estimate of the \eqn{K} function
    or one of its variants.
    An object of class \code{"fv"}.
  }
  \item{\dots}{
    Arguments controlling the smoothing spline
    function \code{smooth.spline}.
  }
  \item{method}{
    Letter \code{"a"}, \code{"b"}, \code{"c"} or \code{"d"} indicating the
    method for deriving the pair correlation function from the
    \code{K} function.
  }
} 
\value{
  A function value table
  (object of class \code{"fv"}, see \code{\link{fv.object}})
  representing a pair correlation function.

  Essentially a data frame containing (at least) the variables
  \item{r}{the vector of values of the argument \eqn{r} 
    at which the pair correlation function \eqn{g(r)} has been  estimated
  }
  \item{pcf}{vector of values of \eqn{g(r)}
  }
}
\details{
  The pair correlation function of a stationary point process is
  \deqn{
    g(r) = \frac{K'(r)}{2\pi r}
  }{
    g(r) = K'(r)/ ( 2 * pi * r) 
  }
  where \eqn{K'(r)} is the derivative of \eqn{K(r)}, the
  reduced second moment function (aka ``Ripley's \eqn{K} function'')
  of the point process. See \code{\link{Kest}} for information
  about \eqn{K(r)}. For a stationary Poisson process, the
  pair correlation function is identically equal to 1. Values
  \eqn{g(r) < 1} suggest inhibition between points;
  values greater than 1 suggest clustering.

  We also apply the same definition to
  other variants of the classical \eqn{K} function,
  such as the multitype \eqn{K} functions
  (see \code{\link{Kcross}}, \code{\link{Kdot}}) and the
  inhomogeneous \eqn{K} function (see \code{\link{Kinhom}}).
  For all these variants, the benchmark value of
  \eqn{K(r) = \pi r^2}{K(r) = pi * r^2} corresponds to
  \eqn{g(r) = 1}.

  This routine computes an estimate of \eqn{g(r)}
  from an estimate of \eqn{K(r)} or its variants,
  using smoothing splines to approximate the derivative.
  It is a method for the generic function \code{\link{pcf}}
  for the class \code{"fv"}.
  
  The argument \code{X} should be an estimated \eqn{K} function,
  given as a function value table (object of class \code{"fv"},
  see \code{\link{fv.object}}).
  This object should be the value returned by
  \code{\link{Kest}}, \code{\link{Kcross}}, \code{\link{Kmulti}}
  or \code{\link{Kinhom}}.
  
  The smoothing spline operations are performed by
  \code{\link{smooth.spline}} and \code{\link{predict.smooth.spline}}
  from the \code{modreg} library.
  Four numerical methods are available:
  \itemize{
    \item
    \bold{"a"} apply smoothing to \eqn{K(r)},
    estimate its derivative, and plug in to the formula above;
    \item 
    \bold{"b"} apply smoothing to
    \eqn{Y(r) = \frac{K(r)}{2 \pi r}}{Y(r) = K(r)/(2 * pi * r)}
    constraining \eqn{Y(0) = 0},
    estimate the derivative of \eqn{Y}, and solve;
    \item
    \bold{"c"} apply smoothing to 
    \eqn{Z(r) = \frac{K(r)}{\pi r^2}}{Y(r) = K(r)/(pi * r^2)}
    constraining \eqn{Z(0)=1},
    estimate its derivative, and solve.
    \item
    \bold{"d"} apply smoothing to 
    \eqn{V(r) = \sqrt{K(r)}}{V(r) = sqrt(K(r))},
    estimate its derivative, and solve.
  }
  Method \code{"c"} seems to be the best at 
  suppressing variability for small values of \eqn{r}.
  However it effectively constrains \eqn{g(0) = 1}.
  If the point pattern seems to have inhibition at small distances,
  you may wish to experiment with method \code{"b"} which effectively
  constrains \eqn{g(0)=0}. Method \code{"a"} seems
  comparatively unreliable.

  Useful arguments to control the splines
  include the smoothing tradeoff parameter \code{spar}
  and the degrees of freedom \code{df}. See \code{\link{smooth.spline}}
  for details.
}
\references{
  Stoyan, D, Kendall, W.S. and Mecke, J. (1995)
  \emph{Stochastic geometry and its applications}.
  2nd edition. Springer Verlag.

  Stoyan, D. and Stoyan, H. (1994)
  Fractals, random shapes and point fields:
  methods of geometrical statistics.
  John Wiley and Sons.
}
\seealso{
  \code{\link{pcf}},
  \code{\link{pcf.ppp}},
  \code{\link{Kest}},
  \code{\link{Kinhom}},
  \code{\link{Kcross}},
  \code{\link{Kdot}},
  \code{\link{Kmulti}},
  \code{\link{alltypes}},
  \code{\link{smooth.spline}},
  \code{\link{predict.smooth.spline}}
}
\examples{
  # univariate point pattern
  data(simdat)
  \testonly{
    simdat <- simdat[seq(1,simdat$n, by=4)]
  }
  K <- Kest(simdat)
  p <- pcf.fv(K, spar=0.5, method="b")
  plot(p, main="pair correlation function for simdat")
  # indicates inhibition at distances r < 0.3
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{nonparametric}
