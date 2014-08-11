\name{Fest}
\alias{Fest}
\alias{empty.space}
\title{Estimate the empty space function F}
\description{
  Estimates the empty space function \eqn{F(r)} from a point pattern in a 
  window of arbitrary shape.
}
\usage{
Fest(X, ..., eps, r=NULL, breaks=NULL, correction=c("rs", "km", "cs"))
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of \eqn{F(r)} will be computed.
    An object of class \code{ppp}, or data
    in any format acceptable to \code{\link{as.ppp}()}.
  }
  \item{\dots}{Ignored.}
  \item{eps}{Optional. A positive number.
    The resolution of the discrete approximation to Euclidean
    distance (see below). There is a sensible default.
  }
  \item{r}{Optional. Numeric vector. The values of the argument \eqn{r}
    at which \eqn{F(r)} should be evaluated. There is a sensible default.
    First-time users are strongly advised not to specify this argument.
    See below for important conditions on \eqn{r}.
  }
  \item{breaks}{
    This argument is for internal use only.
  }
  \item{correction}{
    Optional.
    The edge correction(s) to be used to estimate \eqn{F(r)}.
    A vector of character strings selected from
    \code{"none"}, \code{"rs"}, \code{"km"}, \code{"cs"}
    and \code{"best"}.
  }
}
\value{
  An object of class \code{"fv"}, see \code{\link{fv.object}},
  which can be plotted directly using \code{\link{plot.fv}}.

  Essentially a data frame containing up to seven columns:
  \item{r}{the values of the argument \eqn{r} 
    at which the function \eqn{F(r)} has been  estimated
  }
  \item{rs}{the ``reduced sample'' or ``border correction''
    estimator of \eqn{F(r)}
  }
  \item{km}{the spatial Kaplan-Meier estimator of \eqn{F(r)}
  }
  \item{hazard}{the hazard rate \eqn{\lambda(r)}{lambda(r)}
    of \eqn{F(r)} by the spatial Kaplan-Meier method
  }
  \item{cs}{the Chiu-Stoyan estimator of \eqn{F(r)}
  }
  \item{raw}{the uncorrected estimate of \eqn{F(r)},
  i.e. the empirical distribution of the distance from 
  a random point in the window to the nearest point of
  the data pattern \code{X}
  }
  \item{theo}{the theoretical value of \eqn{F(r)}
  for a stationary Poisson process of the same estimated intensity.
  }
}
\details{
  The empty space function 
  (also called the ``\emph{spherical contact distribution}''
  or the ``\emph{point-to-nearest-event}'' distribution)
  of a stationary point process \eqn{X}
  is the cumulative distribution function \eqn{F} of the distance
  from a fixed point in space to the nearest point of \eqn{X}.

  An estimate of \eqn{F} derived from a spatial point pattern dataset
  can be used in exploratory data analysis and formal inference
  about the pattern (Cressie, 1991; Diggle, 1983; Ripley, 1988).
  In exploratory analyses, the estimate of \eqn{F} is a useful statistic 
  summarising the sizes of gaps in the pattern.
  For inferential purposes, the estimate of \eqn{F} is usually compared to the 
  true value of \eqn{F} for a completely random (Poisson) point process,
  which is
  \deqn{F(r) = 1 - e^{ - \lambda \pi r^2}}{%
    F(r) = 1 - exp( - lambda * pi * r^2)      %
  }
  where \eqn{\lambda}{lambda}
  is the intensity (expected number of points per unit area).
  Deviations between the empirical and theoretical \eqn{F} curves
  may suggest spatial clustering or spatial regularity.

  This algorithm estimates the empty space function \eqn{F} 
  from the point pattern \code{X}. It assumes that \code{X} can be treated
  as a realisation of a stationary (spatially homogeneous) 
  random spatial point process in the plane, observed through
  a bounded window. 
  The window (which is specified in \code{X}) may have arbitrary shape. 

  The argument \code{X} is interpreted as a point pattern object 
  (of class \code{"ppp"}, see \code{\link{ppp.object}}) and can
  be supplied in any of the formats recognised
  by \code{\link{as.ppp}}.

  The algorithm uses two discrete approximations which are controlled
  by the parameter \code{eps} and by the spacing of values of \code{r}
  respectively. (See below for details.)
  First-time users are strongly advised not to specify these arguments.

  The estimation of \eqn{F} is hampered by edge effects arising from 
  the unobservability of points of the random pattern outside the window. 
  An edge correction is needed to reduce bias (Baddeley, 1998; Ripley, 1988). 
  The edge corrections implemented here are the border method or
  "\emph{reduced sample}" estimator, the spatial Kaplan-Meier estimator
  (Baddeley and Gill, 1997) and the Chiu-Stoyan estimator (Chiu and
  Stoyan, 1998). 

  Our implementation makes essential use of the distance transform algorithm
  of image processing (Borgefors, 1986). A fine grid of pixels is 
  created in the observation window. The Euclidean distance between two pixels
  is approximated by the length of the shortest path joining them in the grid,
  where a path is a sequence of steps between adjacent pixels, and 
  horizontal, vertical and diagonal steps have length
  \eqn{1}, \eqn{1} and \eqn{\sqrt 2}{sqrt(2)}
  respectively in pixel units. If the pixel grid is sufficiently fine then
  this is an accurate approximation. 

  The parameter \code{eps}
  is the pixel width of the rectangular raster
  used to compute the distance transform (see below). It must not be too
  large: the absolute error in distance values due to discretisation is bounded
  by \code{eps}.

  If \code{eps} is not specified, the function
  checks whether the window \code{X$window} contains pixel raster
  information. If so, then \code{eps} is set equal to the 
  pixel width of the raster; otherwise, \code{eps}
  defaults to 1/100 of the width of the observation window.

  The argument \code{r} is the vector of values for the
  distance \eqn{r} at which \eqn{F(r)} should be evaluated. 
  It is also used to determine the breakpoints
  (in the sense of \code{\link{hist}}
  for the computation of histograms of distances. The 
  estimators are computed from histogram counts. 
  This introduces a discretisation
  error which is controlled by the fineness of the breakpoints.

  First-time users would be strongly advised not to specify \code{r}.
  However, if it is specified, \code{r} must satisfy \code{r[1] = 0}, 
  and \code{max(r)} must be larger than the radius of the largest disc 
  contained in the window. Furthermore, the spacing of successive
  \code{r} values must be very fine (ideally not greater than \code{eps/4}).

  The algorithm also returns an estimate of the hazard rate function, 
  \eqn{\lambda(r)}{lambda(r)}, of \eqn{F(r)}. The hazard rate is
  defined by
  \deqn{\lambda(r) = - \frac{d}{dr} \log(1 - F(r))}{%
    lambda(r) = - (d/dr) log(1 - F(r)) %
  }
  The hazard rate of \eqn{F} has been proposed as a useful
  exploratory statistic (Baddeley and Gill, 1994).
  The estimate of \eqn{\lambda(r)}{lambda(r)} given here
  is a discrete approximation to the hazard rate of the 
  Kaplan-Meier estimator of \eqn{F}. Note that \eqn{F} is 
  absolutely continuous (for any stationary point process \eqn{X}), 
  so the hazard function always exists (Baddeley and Gill, 1997). 

  The naive empirical distribution of distances from each location
  in the window to the nearest point of the data pattern, is a biased
  estimate of \eqn{F}. However this is also returned by the algorithm
  (if \code{correction="none"}),
  as it is sometimes useful in other contexts.
  Care should be taken not to use the uncorrected
  empirical \eqn{F} as if it were an unbiased estimator of \eqn{F}.
}
\note{
  Sizeable amounts of memory may be needed during the calculation.
}
\references{
  Baddeley, A.J. Spatial sampling and censoring.
     In O.E. Barndorff-Nielsen, W.S. Kendall and
     M.N.M. van Lieshout (eds) 
     \emph{Stochastic Geometry: Likelihood and Computation}.
     Chapman and Hall, 1998.
     Chapter 2, pages 37-78.
  
  Baddeley, A.J. and Gill, R.D. 
    The empty space hazard of a spatial pattern.
    Research Report 1994/3, Department of Mathematics,
    University of Western Australia, May 1994.

  Baddeley, A.J. and Gill, R.D.
     Kaplan-Meier estimators of interpoint distance
		distributions for spatial point processes.
     \emph{Annals of Statistics} \bold{25} (1997) 263-292.

  Borgefors, G.
     Distance transformations in digital images.
     \emph{Computer Vision, Graphics and Image Processing}
     \bold{34} (1986) 344-371.

  Chiu, S.N. and Stoyan, D. (1998)
  Estimators of distance distributions for spatial patterns.
  \emph{Statistica Neerlandica} \bold{52}, 239--246.
  
  Cressie, N.A.C. \emph{Statistics for spatial data}.
    John Wiley and Sons, 1991.

  Diggle, P.J. \emph{Statistical analysis of spatial point patterns}.
    Academic Press, 1983.

  Ripley, B.D. \emph{Statistical inference for spatial processes}.
  Cambridge University Press, 1988.

  Stoyan, D, Kendall, W.S. and Mecke, J.
  \emph{Stochastic geometry and its applications}.
  2nd edition. Springer Verlag, 1995.
}
\section{Warnings}{
  The reduced sample (border method)
  estimator of \eqn{F} is pointwise approximately 
  unbiased, but need not be a valid distribution function; it may 
  not be a nondecreasing function of \eqn{r}. Its range is always 
  within \eqn{[0,1]}.

  The spatial Kaplan-Meier estimator of \eqn{F} is always nondecreasing
  but its maximum value may be less than \eqn{1}.

  The estimate of \eqn{\lambda(r)}{lambda(r)}
  returned by the algorithm is an approximately
  unbiased estimate for the integral of \eqn{\lambda()}{lambda()}
  over the corresponding histogram cell.
  It may exhibit oscillations due to discretisation effects.
  We recommend modest smoothing, such as kernel smoothing with 
  kernel width equal to the width of a histogram cell.
}
\seealso{
  \code{\link{Gest}},
  \code{\link{Jest}},
  \code{\link{Kest}},
  \code{\link{km.rs}},
  \code{\link{reduced.sample}},
  \code{\link{kaplan.meier}}
}
\examples{
   data(cells)
   Fc <- Fest(cells, 0.01)

   # Tip: don't use F for the left hand side!
   # That's an abbreviation for FALSE

   plot(Fc)

   # P-P style plot
   plot(Fc, cbind(km, theo) ~ theo)

   # The empirical F is above the Poisson F
   # indicating an inhibited pattern

   \dontrun{
   plot(Fc, . ~ theo)
   plot(Fc, asin(sqrt(.)) ~ asin(sqrt(theo)))
   }
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{nonparametric}

