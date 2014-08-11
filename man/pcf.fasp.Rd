\name{pcf.fasp}
\alias{pcf.fasp}
\title{Pair Correlation Function obtained from array of K functions}
\description{
  Estimates the (bivariate) pair correlation functions of
  a point pattern, given an array of (bivariate) K functions.
}
\usage{
 \method{pcf}{fasp}(X, \dots, method="c")
}
\arguments{
  \item{X}{
    An array of multitype \eqn{K} functions
    (object of class \code{"fasp"}).
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
  A function array (object of class \code{"fasp"},
  see \code{\link{fasp.object}})
  representing an array of pair correlation functions.
  This can be thought of as a matrix \code{Y} each of whose entries
  \code{Y[i,j]} is a function value table (class \code{"fv"})
  representing the pair correlation function between
  points of type \code{i} and points of type \code{j}.
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
  from an array of estimates of \eqn{K(r)} or its variants,
  using smoothing splines to approximate the derivatives.
  It is a method for the generic function \code{\link{pcf}}.

  The argument \code{X} should be
  a function array (object of class \code{"fasp"},
  see \code{\link{fasp.object}})
  containing several estimates of \eqn{K} functions.
  This should have been obtained from \code{\link{alltypes}}
  with the argument \code{fun="K"}.
  
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
  # multitype point pattern
  KK <- alltypes(amacrine, "K")
  p <- pcf.fasp(KK, spar=0.5, method="b")
  plot(p)
  # strong inhibition between points of the same type
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{nonparametric}
