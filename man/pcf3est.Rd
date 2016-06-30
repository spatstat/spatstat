\name{pcf3est}
\Rdversion{1.1}
\alias{pcf3est}
\title{
  Pair Correlation Function of a Three-Dimensional Point Pattern
}
\description{
  Estimates the pair correlation function
  from a three-dimensional point pattern.
}
\usage{
pcf3est(X, ..., rmax = NULL, nrval = 128, correction = c("translation",
"isotropic"), delta=NULL, adjust=1, biascorrect=TRUE)
}
\arguments{
  \item{X}{
    Three-dimensional point pattern (object of class \code{"pp3"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{rmax}{
    Optional. Maximum value of argument \eqn{r} for which
    \eqn{g_3(r)}{g3(r)} will be estimated. 
  }
  \item{nrval}{
    Optional. Number of values of \eqn{r} for which
    \eqn{g_3(r)}{g3(r)} will be estimated. 
  }
  \item{correction}{
    Optional. Character vector specifying the edge correction(s)
    to be applied. See Details.
  }
  \item{delta}{
    Optional. Half-width of the Epanechnikov smoothing kernel.
  }
  \item{adjust}{
    Optional. Adjustment factor for the default value of \code{delta}.
  }
  \item{biascorrect}{
    Logical value. Whether to correct for underestimation due to
    truncation of the kernel near \eqn{r=0}.
  }
}
\details{
  For a stationary point process \eqn{\Phi}{Phi} in three-dimensional
  space, the pair correlation function is
  \deqn{
    g_3(r) = \frac{K_3'(r)}{4\pi r^2}
  }{
    g3(r) = K3'(r)/(4 * pi * r^2)
  }
  where \eqn{K_3'}{K3'} is the derivative of the
  three-dimensional \eqn{K}-function (see \code{\link{K3est}}).
  
  The three-dimensional point pattern \code{X} is assumed to be a
  partial realisation of a stationary point process \eqn{\Phi}{Phi}.
  The distance between each pair of distinct points is computed.
  Kernel smoothing is applied to these distance values (weighted by
  an edge correction factor) and the result is 
  renormalised to give the estimate of \eqn{g_3(r)}{g3(r)}.

  The available edge corrections are:
  \describe{
    \item{\code{"translation"}:}{
      the Ohser translation correction estimator
      (Ohser, 1983; Baddeley et al, 1993)
    }
    \item{\code{"isotropic"}:}{
      the three-dimensional counterpart of
      Ripley's isotropic edge correction (Ripley, 1977; Baddeley et al, 1993).
    }
  }

  Kernel smoothing is performed using the Epanechnikov kernel
  with half-width \code{delta}. If \code{delta} is missing, the
  default is to use the rule-of-thumb
  \eqn{\delta = 0.26/\lambda^{1/3}}{delta = 0.26/lambda^(1/3)} where
  \eqn{\lambda = n/v}{lambda = n/v} is the estimated intensity, computed
  from the number \eqn{n} of data points and the volume \eqn{v} of the
  enclosing box. This default value of \code{delta} is multiplied by
  the factor \code{adjust}. 

  The smoothing estimate of the pair correlation \eqn{g_3(r)}{g3(r)}
  is typically an underestimate when \eqn{r} is small, due to
  truncation of the kernel at \eqn{r=0}. 
  If \code{biascorrect=TRUE}, the smoothed estimate is
  approximately adjusted for this bias. This is advisable whenever
  the dataset contains a sufficiently large number of points.
}
\value{
  A function value table (object of class \code{"fv"}) that can be
  plotted, printed or coerced to a data frame containing the function
  values.

  Additionally the value of \code{delta} is returned as an attribute
  of this object.
}
\references{
  Baddeley, A.J, Moyeed, R.A., Howard, C.V. and Boyde, A. (1993)
  Analysis of a three-dimensional point pattern with replication.
  \emph{Applied Statistics} \bold{42}, 641--668.

  Ohser, J. (1983)
  On estimators for the reduced second moment measure of
  point processes. \emph{Mathematische Operationsforschung und
  Statistik, series Statistics}, \bold{14}, 63 -- 71.

  Ripley, B.D. (1977)
  Modelling spatial patterns (with discussion).
  \emph{Journal of the Royal Statistical Society, Series B},
  \bold{39}, 172 -- 212.
}
\author{
  \adrian
  
  
  and Rana Moyeed.
}
\seealso{
  \code{\link{K3est}},
  \code{\link{pcf}}
}
\examples{
  X <- rpoispp3(250)
  Z <- pcf3est(X)
  Zbias <- pcf3est(X, biascorrect=FALSE)
  if(interactive()) {
    opa <- par(mfrow=c(1,2))
    plot(Z,     ylim.covers=c(0, 1.2))
    plot(Zbias, ylim.covers=c(0, 1.2))
    par(opa)
  }
  attr(Z, "delta")
}
\keyword{spatial}
\keyword{nonparametric}
