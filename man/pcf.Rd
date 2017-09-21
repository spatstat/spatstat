\name{pcf}
\alias{pcf}
\title{Pair Correlation Function}
\description{
  Estimate the pair correlation function.
}
\usage{
 pcf(X, \dots)
}
\arguments{
  \item{X}{
    Either the observed data point pattern,
    or an estimate of its \eqn{K} function,
    or an array of multitype \eqn{K} functions
    (see Details).
  }
  \item{\dots}{
    Other arguments passed to the appropriate method.
  }
} 
\value{
  Either a function value table
  (object of class \code{"fv"}, see \code{\link{fv.object}})
  representing a pair correlation function,
  or a function array (object of class \code{"fasp"},
  see \code{\link{fasp.object}})
  representing an array of pair correlation functions.
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
  either directly from a point pattern,
  or indirectly from an estimate of \eqn{K(r)} or one of its variants.

  This function is generic, with methods for
  the classes \code{"ppp"}, \code{"fv"} and \code{"fasp"}.

  If \code{X} is a point pattern (object of class \code{"ppp"})
  then the pair correlation function is estimated using
  a traditional kernel smoothing method (Stoyan and Stoyan, 1994).
  See \code{\link{pcf.ppp}} for details.

  If \code{X} is a function value table (object of class \code{"fv"}),
  then it is assumed to contain estimates of the \eqn{K} function
  or one of its variants (typically obtained from \code{\link{Kest}} or
  \code{\link{Kinhom}}).
  This routine computes an estimate of \eqn{g(r)} 
  using smoothing splines to approximate the derivative.
  See \code{\link{pcf.fv}} for details.

  If \code{X} is a function value array (object of class \code{"fasp"}),
  then it is assumed to contain estimates of several \eqn{K} functions
  (typically obtained from \code{\link{Kmulti}} or
  \code{\link{alltypes}}). This routine computes
  an estimate of \eqn{g(r)} for each cell in the array,
  using smoothing splines to approximate the derivatives.
  See \code{\link{pcf.fasp}} for details.
}
\references{
  Stoyan, D. and Stoyan, H. (1994)
  Fractals, random shapes and point fields:
  methods of geometrical statistics.
  John Wiley and Sons.
}
\seealso{
  \code{\link{pcf.ppp}},
  \code{\link{pcf.fv}},
  \code{\link{pcf.fasp}},
  \code{\link{Kest}},
  \code{\link{Kinhom}},
  \code{\link{Kcross}},
  \code{\link{Kdot}},
  \code{\link{Kmulti}},
  \code{\link{alltypes}}
}
\examples{
  # ppp object
  X <- simdat
  \testonly{
    X <- X[seq(1,npoints(X), by=4)]
  }
  p <- pcf(X)
  plot(p)

  # fv object
  K <- Kest(X)
  p2 <- pcf(K, spar=0.8, method="b")
  plot(p2)

  # multitype pattern; fasp object
  amaK <- alltypes(amacrine, "K")
  amap <- pcf(amaK, spar=1, method="b")
  plot(amap)
}
\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{nonparametric}
