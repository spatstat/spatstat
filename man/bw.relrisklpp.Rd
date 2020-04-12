\name{bw.relrisklpp}
\alias{bw.relrisklpp}
\title{
  Cross Validated Bandwidth Selection for Relative Risk Estimation
  on a Network
}
\description{
  Uses cross-validation to select a smoothing bandwidth
  for the estimation of relative risk on a linear network.
}
\usage{
  bw.relrisklpp(X, \dots,
    method = c("likelihood", "leastsquares", "KelsallDiggle", "McSwiggan"),
    distance=c("path", "euclidean"),
    hmin = NULL, hmax = NULL, nh = NULL,
    fast = TRUE, fastmethod = "onestep",
    floored = TRUE, reference = c("thumb", "uniform", "sigma"),
    allow.infinite = TRUE, epsilon = 1e-20, fudge = 0,
    verbose = FALSE, warn = TRUE)
}
\arguments{
  \item{X}{
    A multitype point pattern on a linear network (object of class
    \code{"lpp"} which has factor-valued marks).
  }
  \item{\dots}{
    Arguments passed to \code{\link{density.lpp}} to control the
    resolution of the algorithm.
  }
  \item{method}{
    Character string (partially matched)
    determining the cross-validation method.
    See Details.
  }
  \item{distance}{
    Character string (partially matched)
    specifying the type of smoothing kernel.
    See \code{\link{density.lpp}}.
  }
  \item{hmin,hmax}{
    Optional. Numeric values.
    Range of trial values of smoothing bandwith \code{sigma}
    to consider. There is a sensible default.
  }
  \item{nh}{
    Number of trial values of smoothing bandwidth \code{sigma}
    to consider.
  }
  \item{fast}{
    Logical value specifying whether the leave-one-out density estimates
    should be computed using a fast approximation (\code{fast=TRUE}, the
    default) or exactly (\code{fast=FALSE}).
  }
  \item{fastmethod, floored}{
    Developer use only.
  }
  \item{reference}{
    Character string (partially matched) specifying the
    bandwidth for calculating the
    reference intensities used in the McSwiggan method
    (modified Kelsall-Diggle method).
    \code{reference="sigma"} means the maximum bandwidth considered,
    which is given by the argument \code{sigma}.
    \code{reference="thumb"} means the bandwidths selected by
    Scott's rule of thumb \code{\link{bw.scott.iso}}.
    \code{reference="uniform"} means infinite bandwidth corresponding to
    uniform intensity.
  }
  \item{allow.infinite}{
    Logical value indicating whether an infinite bandwidth
    (corresponding to a constant relative risk) should be permitted
    as a possible choice of bandwidth.
  }
  \item{epsilon}{
    A small constant value added to the reference density in some
    of the cross-validation calculations, to improve performance.
  }
  \item{fudge}{
    Fudge factor to prevent very small density estimates in the
    leave-one-out calculation. If \code{fudge > 0},
    then the lowest permitted value
    for a leave-one-out estimate of intensity is
    \code{fudge/L}, where \code{L} is the total length of the
    network. 
  }
  \item{verbose}{
    Logical value indicating whether to print progress reports,
  }
  \item{warn}{
    Logical. If \code{TRUE}, issue a warning if the minimum of
    the cross-validation criterion occurs at one of the ends of the
    search interval.
  }
}
\details{
  This function computes an optimal value of smoothing bandwidth
  for the nonparametric estimation of relative risk on a linear network
  using \code{\link{relrisk.lpp}}.
  The optimal value is found by optimising a cross-validation criterion.

  The cross-validation criterion is selected by the argument \code{method}:
  \tabular{ll}{
    \code{method="likelihood"} \tab likelihood cross-validation \cr
    \code{method="leastsquares"} \tab least squares cross-validation \cr
    \code{method="KelsallDiggle"} \tab Kelsall and Diggle (1995) density
  ratio cross-validation \cr
    \code{method="McSwiggan"} \tab McSwiggan et al (2019) modified
  density ratio cross-validation \cr
  }
  See McSwiggan et al (2019) for details.

  The result is a numerical value giving the selected bandwidth \code{sigma}.
  The result also belongs to the class \code{"bw.optim"}
  allowing it to be printed and plotted. The plot shows the cross-validation
  criterion as a function of bandwidth.  

  The range of values for the smoothing bandwidth \code{sigma}
  is set by the arguments \code{hmin, hmax}. There is a sensible default,
  based on the linear network version of Scott's rule
  \code{\link{bw.scott.iso}}.
  
  If the optimal bandwidth is achieved at an endpoint of the
  interval \code{[hmin, hmax]}, the algorithm will issue a warning
  (unless \code{warn=FALSE}). If this occurs, then it is probably advisable
  to expand the interval by changing the arguments \code{hmin, hmax}.

  The cross-validation procedure is based on kernel estimates
  of intensity, which are computed by \code{\link{density.lpp}}.
  Any arguments \code{...} are passed to \code{\link{density.lpp}}
  to control the kernel estimation procedure. This includes the
  argument \code{distance} which specifies the type of kernel.
  The default is \code{distance="path"};
  the fastest option is \code{distance="euclidean"}.
}
\value{
  A numerical value giving the selected bandwidth.
  The result also belongs to the class \code{"bw.optim"}
  which can be plotted.
}
\references{
  Kelsall, J.E. and Diggle, P.J. (1995)
  Kernel estimation of relative risk.
  \emph{Bernoulli} \bold{1}, 3--16.

  McSwiggan, G., Baddeley, A. and Nair, G. (2019)
  Estimation of relative risk for events on a linear network.
  \emph{Statistics and Computing} \bold{30} (2) 469--484.
}
\author{
  Greg McSwiggan and \adrian.
}
\seealso{
  \code{\link{relrisk.lpp}}
}
\examples{
   set.seed(2020)
   X <- superimpose(A=runiflpp(20, simplenet),
                    B=runifpointOnLines(20, as.psp(simplenet)[1]))
   plot(bw.relrisklpp(X, hmin=0.1, hmax=0.3, method="McSwiggan"))
   plot(bw.relrisklpp(X, hmin=0.1, hmax=0.3, distance="euclidean"))
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}

