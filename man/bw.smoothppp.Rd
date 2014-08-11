\name{bw.smoothppp}
\alias{bw.smoothppp}
\title{
  Cross Validated Bandwidth Selection for Spatial Smoothing
}
\description{
  Uses least-squares cross-validation to select a smoothing bandwidth
  for spatial smoothing of marks.
}
\usage{
   bw.smoothppp(X, nh = spatstat.options("n.bandwidth"),
   hmin=NULL, hmax=NULL, warn=TRUE)
}

\arguments{
  \item{X}{
    A marked point pattern with numeric marks.
  }
  \item{nh}{
    Number of trial values of smoothing bandwith \code{sigma}
    to consider. The default is 32.
  }
  \item{hmin, hmax}{
    Optional. Numeric values.
    Range of trial values of smoothing bandwith \code{sigma}
    to consider. There is a sensible default.
  }
  \item{warn}{
    Logical. If \code{TRUE}, issue a warning if the minimum of
    the cross-validation criterion occurs at one of the ends of the
    search interval.
  }
}
\details{
  This function selects an appropriate bandwidth for the nonparametric
  smoothing of mark values using \code{\link{Smooth.ppp}}.
  
  The argument \code{X} must be a marked point pattern
  with a vector or data frame of marks. All mark values must be numeric.
  
  The bandwidth is selected by least-squares cross-validation.
  Let \eqn{y_i}{y[i]} be the mark value at the \eqn{i}th data point.
  For a particular choice of smoothing bandwidth,
  let \eqn{\hat y_i}{y*[i]} be the smoothed value at the \eqn{i}th data point.
  Then the bandwidth is chosen to minimise 
  the squared error of the smoothed values
  \eqn{\sum_i (y_i - \hat y_i)^2}{sum (y[i] - y*[i])^2}.

  The result of \code{bw.smoothppp}
  is a numerical value giving the selected bandwidth \code{sigma}.
  The result also belongs to the class \code{"bw.optim"}
  allowing it to be printed and plotted. The plot shows the cross-validation
  criterion as a function of bandwidth.  
  
  The range of values for the smoothing bandwidth \code{sigma}
  is set by the arguments \code{hmin, hmax}. There is a sensible default,
  based on the nearest neighbour distances.
  
  If the optimal bandwidth is achieved at an endpoint of the
  interval \code{[hmin, hmax]}, the algorithm will issue a warning
  (unless \code{warn=FALSE}). If this occurs, then it is probably advisable
  to expand the interval by changing the arguments \code{hmin, hmax}.

  Computation time depends on the number \code{nh} of trial values
  considered, and also on the range \code{[hmin, hmax]} of values
  considered, because larger values of \code{sigma} require
  calculations involving more pairs of data points.
}
\value{
  A numerical value giving the selected bandwidth.
  The result also belongs to the class \code{"bw.optim"}
  which can be plotted.
}
\seealso{
  \code{\link{Smooth.ppp}}
}
\examples{
  data(longleaf)
  \testonly{op <- spatstat.options(n.bandwidth=8)}
  b <- bw.smoothppp(longleaf)
  b
  plot(b)
  \testonly{spatstat.options(op)}
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
