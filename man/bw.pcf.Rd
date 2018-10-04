\name{bw.pcf}
\alias{bw.pcf}
\title{
  Cross Validated Bandwidth Selection for Pair Correlation Function
}
\description{
  Uses composite likelihood or generalized least squares 
  cross-validation to select a smoothing bandwidth
  for the kernel estimation of pair correlation function.
}
\usage{
  bw.pcf(X, rmax=NULL, lambda=NULL, divisor="r", 
         kernel="epanechnikov", nr=10000, bias.correct=TRUE, 
         cv.method=c("compLik", "leastSQ"), simple=TRUE, srange=NULL,
	 \dots, verbose=FALSE)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
  \item{rmax}{
    Numeric. Maximum value of the spatial lag distance \eqn{r} 
    for which \eqn{g(r)} should be evaluated.
  }
  \item{lambda}{
    Optional.
    Values of the estimated intensity function.
    A vector giving the intensity values
    at the points of the pattern \code{X}.
  }
  \item{divisor}{
    Choice of divisor in the estimation formula:
    either \code{"r"} (the default) or \code{"d"}. 
    See \code{pcf.ppp}.
  }
  \item{kernel}{
    Choice of smoothing kernel, passed to \code{density}; 
    see \code{\link{pcf}} and \code{\link{pcfinhom}}.
  }
  \item{nr}{
    Integer. Number of subintervals for discretization of 
    [0, rmax] to use in computing numerical integrals.
  }
  \item{bias.correct}{
    Logical. Whether to use bias corrected version of the kernel 
    estimate. See Details.
  }
  \item{cv.method}{
    Choice of cross validation method: either
    \code{"compLik"} or \code{"leastSQ"} (partially matched).
  }
  \item{simple}{
    Logical. Whether to use simple removal of spatial lag 
    distances. See Details.
  }
  \item{srange}{
    Optional. Numeric vector of length 2 giving the range of
    bandwidth values that should be searched to find the optimum
    bandwidth.
  }
  \item{\dots}{
    Other arguments, passed to \code{\link{pcf}} or 
    \code{\link{pcfinhom}}.
  }
  \item{verbose}{
    Logical value indicating whether to print progress reports
    during the optimization procedure.
  }
}
\details{
  This function selects an appropriate bandwidth \code{bw}
  for the kernel estimator of the pair correlation function 
  of a point process intensity computed by \code{\link{pcf.ppp}} 
  (homogeneous case) or \code{\link{pcfinhom}} 
  (inhomogeneous case).

  With \code{cv.method="leastSQ"}, the bandwidth 
  \eqn{h} is chosen to minimise an unbiased 
  estimate of the integrated mean-square error criterion 
  \eqn{M(h)} defined in equation (4) in Guan (2007a).
  The code implements the fast algorithm of Jalilian and Waagepetersen
  (2018).
    
  With \code{cv.method="compLik"}, the bandwidth 
  \eqn{h} is chosen to maximise a likelihood 
  cross-validation criterion \eqn{CV(h)} defined in 
  equation (6) of Guan (2007b).
  
  \deqn{
    M(b) = \frac{\mbox{MSE}(\sigma)}{\lambda^2} - g(0)
  }{
    M(b) = \int_{0}^{rmax} \hat{g}^2(r;b) r dr - \sum_{u,v}
  }

  The result is a numerical value giving the selected bandwidth.
}
\section{Definition of bandwidth}{
  The bandwidth \code{bw} returned by \code{bw.pcf}
  corresponds to the standard deviation of the smoothoing 
  kernel. As mentioned in the documentation of 
  \code{\link{density.default}} and \code{\link{pcf.ppp}}, 
  this differs from the scale parameter \code{h} of 
  the smoothing kernel which is often considered in the 
  literature as the bandwidth of the kernel function.
  For example for the Epanechnikov kernel, \code{bw=h/sqrt(h)}.
}
\value{
  A numerical value giving the selected bandwidth.
  The result also belongs to the class \code{"bw.optim"}
  which can be plotted.
}
\seealso{
  \code{\link{pcf.ppp}},
  \code{\link{pcfinhom}}
}
\examples{
  b <- bw.pcf(redwood)
  plot(pcf(redwood, bw=b))
}
\references{
  Guan, Y. (2007a). 
  A composite likelihood cross-validation approach in selecting 
  bandwidth for the estimation of the pair correlation function. 
  \emph{Scandinavian Journal of Statistics}, 
  \bold{34}(2), 336--346.
  
  Guan, Y. (2007b). 
  A least-squares cross-validation bandwidth selection approach 
  in pair correlation function estimations. 
  \emph{Statistics & Probability Letters}, 
  \bold{77}(18), 1722--1729.

  Jalilian, A. and Waagepetersen, R. (2018)
  Fast bandwidth selection for estimation of the pair correlation
  function.
  \emph{Journal of Statistical Computation and Simulation},
  \bold{88}(10), 2001--2011.
  \url{https://www.tandfonline.com/doi/full/10.1080/00949655.2018.1428606}
}
\author{
  Rasmus Waagepetersen and Abdollah Jalilian. 
  Adapted for \pkg{spatstat} by \spatstatAuthors.
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
