\name{Kest}
\alias{Kest}
\title{K-function}
\description{
Estimates Ripley's reduced second moment function \eqn{K(r)} 
from a point pattern in a window of arbitrary shape.
}
\usage{
  Kest(X, \dots, r=NULL, rmax=NULL, breaks=NULL, 
     correction=c("border", "isotropic", "Ripley", "translate"),
    nlarge=3000, domain=NULL, var.approx=FALSE, ratio=FALSE)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of \eqn{K(r)} will be computed.
    An object of class \code{"ppp"}, or data
    in any format acceptable to \code{\link{as.ppp}()}.
  }
  \item{\dots}{Ignored.}
  \item{r}{
    Optional. Vector of values for the argument \eqn{r} at which \eqn{K(r)} 
    should be evaluated. Users are advised \emph{not} to specify this
    argument; there is a sensible default. If necessary, specify \code{rmax}.
  }
  \item{rmax}{
    Optional. Maximum desired value of the argument \eqn{r}.
  }
  \item{breaks}{
    This argument is for internal use only.
  }
  \item{correction}{
    Optional. A character vector containing any selection of the
    options \code{"none"}, \code{"border"}, \code{"bord.modif"},
    \code{"isotropic"}, \code{"Ripley"}, \code{"translate"},
    \code{"translation"}, \code{"rigid"},
    \code{"none"}, \code{"good"} or \code{"best"}.
    It specifies the edge correction(s) to be applied.
    Alternatively \code{correction="all"} selects all options.
  }
  \item{nlarge}{
    Optional. Efficiency threshold.
    If the number of points exceeds \code{nlarge}, then only the
    border correction will be computed (by default), using a fast algorithm.
  }
  \item{domain}{
    Optional. Calculations will be restricted to this subset
    of the window. See Details.
  }
  \item{var.approx}{Logical. If \code{TRUE}, the approximate
    variance of \eqn{\hat K(r)}{Kest(r)} under CSR
    will also be computed.
  }
  \item{ratio}{
    Logical. 
    If \code{TRUE}, the numerator and denominator of
    each edge-corrected estimate will also be saved,
    for use in analysing replicated point patterns.
  }
}
\value{
  An object of class \code{"fv"}, see \code{\link{fv.object}},
  which can be plotted directly using \code{\link{plot.fv}}.

  Essentially a data frame containing columns
  \item{r}{the vector of values of the argument \eqn{r} 
    at which the function \eqn{K} has been  estimated
  }
  \item{theo}{the theoretical value \eqn{K(r) = \pi r^2}{K(r) = pi * r^2}
    for a stationary Poisson process
  }
  together with columns named 
  \code{"border"}, \code{"bord.modif"},
  \code{"iso"} and/or \code{"trans"},
  according to the selected edge corrections. These columns contain
  estimates of the function \eqn{K(r)} obtained by the edge corrections
  named.

  If \code{var.approx=TRUE} then the return value
  also has columns \code{rip} and \code{ls} containing approximations
  to the variance of \eqn{\hat K(r)}{Kest(r)} under CSR.

  If \code{ratio=TRUE} then the return value also has two
  attributes called \code{"numerator"} and \code{"denominator"}
  which are \code{"fv"} objects
  containing the numerators and denominators of each
  estimate of \eqn{K(r)}. 
}
\details{
  The \eqn{K} function (variously called ``Ripley's K-function''
  and the ``reduced second moment function'')
  of a stationary point process \eqn{X} is defined so that
  \eqn{\lambda K(r)}{lambda K(r)} equals the expected number of
  additional random points within a distance \eqn{r} of a
  typical random point of \eqn{X}. Here \eqn{\lambda}{lambda}
  is the intensity of the process,
  i.e. the expected number of points of \eqn{X} per unit area.
  The \eqn{K} function is determined by the 
  second order moment properties of \eqn{X}.
 
  An estimate of \eqn{K} derived from a spatial point pattern dataset
  can be used in exploratory data analysis and formal inference
  about the pattern (Cressie, 1991; Diggle, 1983; Ripley, 1977, 1988).
  In exploratory analyses, the estimate of \eqn{K} is a useful statistic 
  summarising aspects of inter-point ``dependence'' and ``clustering''.
  For inferential purposes, the estimate of \eqn{K} is usually compared to the 
  true value of \eqn{K} for a completely random (Poisson) point process,
  which is \eqn{K(r) = \pi r^2}{K(r) = pi * r^2}.
  Deviations between the empirical and theoretical \eqn{K} curves
  may suggest spatial clustering or spatial regularity.
 
  This routine \code{Kest} estimates the \eqn{K} function
  of a stationary point process, given observation of the process
  inside a known, bounded window. 
  The argument \code{X} is interpreted as a point pattern object 
  (of class \code{"ppp"}, see \code{\link{ppp.object}}) and can
  be supplied in any of the formats recognised by
  \code{\link{as.ppp}()}.

  The estimation of \eqn{K} is hampered by edge effects arising from 
  the unobservability of points of the random pattern outside the window. 
  An edge correction is needed to reduce bias (Baddeley, 1998; Ripley, 1988). 
  The corrections implemented here are
  \describe{
    \item{border}{the border method or
      ``reduced sample'' estimator (see Ripley, 1988). This is
      the least efficient (statistically) and the fastest to compute.
      It can be computed for a window of arbitrary shape.
    }
    \item{isotropic/Ripley}{Ripley's isotropic correction
      (see Ripley, 1988; Ohser, 1983).
      This is implemented for rectangular and polygonal windows
      (not for binary masks).
    }
    \item{translate/translation}{Translation correction (Ohser, 1983).
      Implemented for all window geometries, but slow for
      complex windows. 
    }
    \item{rigid}{Rigid motion correction (Ohser and Stoyan, 1981).
      Implemented for all window geometries, but slow for
      complex windows. 
    }
    \item{none}{
      Uncorrected estimate.
      An estimate of the K function \emph{without} edge correction.
      (i.e. setting \eqn{e_{ij} = 1}{e[i,j] = 1} in the equation below.
      This estimate is \bold{biased} and should
      not be used for data analysis, \emph{unless} you have
      an extremely large point pattern (more than 100,000 points).
    }
    \item{best}{
      Selects the best edge correction
      that is available for the geometry of the window. Currently
      this is Ripley's isotropic correction for a rectangular
      or polygonal window, and the translation correction for masks.
    }
    \item{good}{
      Selects the best edge correction
      that can be computed in a reasonable time.
      This is the same as \code{"best"} for datasets with fewer
      than 3000 points; otherwise the selected edge correction
      is \code{"border"}, unless there are more than 100,000 points, when 
      it is \code{"none"}.
    }
  }
  
  The estimates of \eqn{K(r)} are of the form
  \deqn{
    \hat K(r) = \frac a {n (n-1) } \sum_i \sum_j I(d_{ij}\le r) e_{ij}
  }{
    Kest(r) = (a/(n * (n-1))) * sum[i,j] I(d[i,j] <= r) e[i,j])
  }
  where \eqn{a} is the area of the window, \eqn{n} is the number of
  data points, and the sum is taken over all ordered pairs of points
  \eqn{i} and \eqn{j} in \code{X}.
  Here \eqn{d_{ij}}{d[i,j]} is the distance between the two points,
  and \eqn{I(d_{ij} \le r)}{I(d[i,j] <= r)} is the indicator
  that equals 1 if the distance is less than or equal to \eqn{r}.
  The term \eqn{e_{ij}}{e[i,j]} is the edge correction weight (which
  depends on the choice of edge correction listed above).

  Note that this estimator assumes the process is stationary (spatially
  homogeneous). For inhomogeneous point patterns, see
  \code{\link{Kinhom}}.

  If the point pattern \code{X} contains more than about 3000 points,
  the isotropic and translation edge corrections can be computationally
  prohibitive. The computations for the border method are much faster,
  and are statistically efficient when there are large numbers of
  points. Accordingly, if the number of points in \code{X} exceeds
  the threshold \code{nlarge}, then only the border correction will be
  computed. Setting \code{nlarge=Inf} or \code{correction="best"}
  will prevent this from happening.
  Setting \code{nlarge=0} is equivalent to selecting only the border
  correction with \code{correction="border"}.

  If \code{X} contains more than about 100,000 points,
  even the border correction is time-consuming. You may want to consider
  setting \code{correction="none"} in this case.
  There is an even faster algorithm for the uncorrected estimate.

  Approximations to the variance of \eqn{\hat K(r)}{Kest(r)}
  are available, for the case of the isotropic edge correction estimator,
  \bold{assuming complete spatial randomness}
  (Ripley, 1988; Lotwick and Silverman, 1982; Diggle, 2003, pp 51-53).
  If \code{var.approx=TRUE}, then the result of
  \code{Kest} also has a column named \code{rip} 
  giving values of Ripley's (1988) approximation to
  \eqn{\mbox{var}(\hat K(r))}{var(Kest(r))},
  and (if the window is a rectangle) a column named \code{ls} giving
  values of Lotwick and Silverman's (1982) approximation.
  
  If the argument \code{domain} is given, the calculations will
  be restricted to a subset of the data. In the formula for \eqn{K(r)} above,
  the \emph{first} point \eqn{i} will be restricted to lie inside
  \code{domain}. The result is an approximately unbiased estimate
  of \eqn{K(r)} based on pairs of points in which the first point lies
  inside \code{domain} and the second point is unrestricted.
  This is useful in bootstrap techniques. The argument \code{domain}
  should be a window (object of class \code{"owin"}) or something acceptable to
  \code{\link{as.owin}}. It must be a subset of the
  window of the point pattern \code{X}.

  The estimator \code{Kest} ignores marks.
  Its counterparts for multitype point patterns
  are \code{\link{Kcross}}, \code{\link{Kdot}},
  and for general marked point patterns
  see \code{\link{Kmulti}}. 

  Some writers, particularly Stoyan (1994, 1995) advocate the use of
  the ``pair correlation function''
  \deqn{
    g(r) = \frac{K'(r)}{2\pi r}
  }{
    g(r) = K'(r)/ ( 2 * pi * r) 
  }
  where \eqn{K'(r)} is the derivative of \eqn{K(r)}.
  See \code{\link{pcf}} on how to estimate this function.
}
\section{Envelopes, significance bands and confidence intervals}{
  To compute simulation envelopes for the \eqn{K}-function
  under CSR, use \code{\link{envelope}}. 
  
  To compute a confidence interval for the true \eqn{K}-function,
  use \code{\link{varblock}} or \code{\link{lohboot}}.
}
\references{
Baddeley, A.J. Spatial sampling and censoring.
     In O.E. Barndorff-Nielsen, W.S. Kendall and
     M.N.M. van Lieshout (eds) 
     \emph{Stochastic Geometry: Likelihood and Computation}.
     Chapman and Hall, 1998.
     Chapter 2, pages 37--78.
  
  Cressie, N.A.C. \emph{Statistics for spatial data}.
    John Wiley and Sons, 1991.

  Diggle, P.J. \emph{Statistical analysis of spatial point patterns}.
  Academic Press, 1983.

  Ohser, J. (1983)
  On estimators for the reduced second moment measure of
  point processes. \emph{Mathematische Operationsforschung und
  Statistik, series Statistics}, \bold{14}, 63 -- 71.
  
  Ohser, J. and Stoyan, D. (1981)
  On the second-order and orientation analysis of
  planar stationary point processes.
  \emph{Biometrical Journal} \bold{23}, 523--533.

  Ripley, B.D. (1977)
  Modelling spatial patterns (with discussion).
  \emph{Journal of the Royal Statistical Society, Series B},
  \bold{39}, 172 -- 212.

  Ripley, B.D. \emph{Statistical inference for spatial processes}.
  Cambridge University Press, 1988.

  Stoyan, D, Kendall, W.S. and Mecke, J. (1995)
  \emph{Stochastic geometry and its applications}.
  2nd edition. Springer Verlag.

  Stoyan, D. and Stoyan, H. (1994)
  Fractals, random shapes and point fields:
  methods of geometrical statistics.
  John Wiley and Sons.
} 
\section{Warnings}{
  The estimator of \eqn{K(r)} is approximately unbiased for each fixed
  \eqn{r}, for point processes which do not have very strong
  interaction. (For point processes with a strong clustering interaction,
  the estimator is negatively biased; for point processes with a strong
  inhibitive interaction, the estimator is positively biased.)

  Bias increases with \eqn{r} and depends on the window geometry.
  For a rectangular window it is prudent to restrict the \eqn{r} values to
  a maximum of \eqn{1/4} of the smaller side length of the rectangle
  (Ripley, 1977, 1988; Diggle, 1983).
  Bias may become appreciable for point patterns consisting of 
  fewer than 15 points.
 
  While \eqn{K(r)} is always a non-decreasing function, the estimator 
  of \eqn{K} is not guaranteed to be non-decreasing. This is rarely 
  a problem in practice, except for the border correction estimators
  when the number of points is small.
}
\seealso{
  \code{\link{localK}} to extract individual summands in the \eqn{K}
  function.

  \code{\link{pcf}} for the pair correlation.

  \code{\link{Fest}},
  \code{\link{Gest}},
  \code{\link{Jest}}
  for alternative summary functions.
  
  \code{\link{Kcross}},
  \code{\link{Kdot}},
  \code{\link{Kinhom}},
  \code{\link{Kmulti}} for counterparts of the \eqn{K} function
  for multitype point patterns.
  
  \code{\link{reduced.sample}} for the calculation of reduced sample
  estimators.
}
\examples{
 X <- runifpoint(50)
 K <- Kest(X)
 K <- Kest(cells, correction="isotropic")
 plot(K)
 plot(K, main="K function for cells")
 # plot the L function
 plot(K, sqrt(iso/pi) ~ r)
 plot(K, sqrt(./pi) ~ r, ylab="L(r)", main="L function for cells")
}
\author{\adrian
  and \rolf
}
\keyword{spatial}
\keyword{nonparametric}

 
 
