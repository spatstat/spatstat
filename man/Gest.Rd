\name{Gest}
\alias{Gest}
\alias{nearest.neighbour}
\title{
  Nearest Neighbour Distance Function G
}
\description{
Estimates the nearest neighbour distance distribution
function \eqn{G(r)} from a point pattern in a 
window of arbitrary shape.
}
\usage{
Gest(X, r=NULL, breaks=NULL, \dots, correction=c("rs", "km", "han"))
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of \eqn{G(r)} will be computed.
    An object of class \code{ppp}, or data
    in any format acceptable to \code{\link{as.ppp}()}.
  }
  \item{r}{Optional. Numeric vector. The values of the argument \eqn{r}
    at which \eqn{G(r)} should be evaluated. There is a sensible default.
    First-time users are strongly advised not to specify this argument.
    See below for important conditions on \eqn{r}.
  }
  \item{breaks}{An alternative to the argument \code{r}.
    Not normally invoked by the user. See the \bold{Details} section.
  }
  \item{\dots}{Ignored.}
  \item{correction}{
    Optional.
    The edge correction(s) to be used to estimate \eqn{G(r)}.
    A vector of character strings selected from
    \code{"none"}, \code{"rs"}, \code{"km"}, \code{"Hanisch"}
    and \code{"best"}.
  }
}
\value{
  An object of class \code{"fv"}, see \code{\link{fv.object}},
  which can be plotted directly using \code{\link{plot.fv}}.

  Essentially a data frame containing some or all of the following
  columns:
  
  \item{r}{the values of the argument \eqn{r} 
    at which the function \eqn{G(r)} has been  estimated
  }
  \item{rs}{the ``reduced sample'' or ``border correction''
    estimator of \eqn{G(r)}
  }
  \item{km}{the spatial Kaplan-Meier estimator of \eqn{G(r)}
  }
  \item{hazard}{the hazard rate \eqn{\lambda(r)}{lambda(r)}
    of \eqn{G(r)} by the spatial Kaplan-Meier method
  }
  \item{raw}{the uncorrected estimate of \eqn{G(r)},
  i.e. the empirical distribution of the distances from 
  each point in the pattern \code{X} to the nearest other point of
  the pattern
  }
  \item{han}{the Hanisch correction estimator of \eqn{G(r)}
  }
  \item{theo}{the theoretical value of \eqn{G(r)}
  for a stationary Poisson process of the same estimated intensity.
  }
}
\details{
  The nearest neighbour distance distribution function 
  (also called the ``\emph{event-to-event}'' or
  ``\emph{inter-event}'' distribution)
  of a point process \eqn{X}
  is the cumulative distribution function \eqn{G} of the distance
  from a typical random point of \eqn{X} to
  the nearest other point of \eqn{X}.

  An estimate of \eqn{G} derived from a spatial point pattern dataset
  can be used in exploratory data analysis and formal inference
  about the pattern (Cressie, 1991; Diggle, 1983; Ripley, 1988).
  In exploratory analyses, the estimate of \eqn{G} is a useful statistic 
  summarising one aspect of the ``clustering'' of points.
  For inferential purposes, the estimate of \eqn{G} is usually compared to the 
  true value of \eqn{G} for a completely random (Poisson) point process,
  which is
  \deqn{G(r) = 1 - e^{ - \lambda \pi r^2} }{%
    G(r) = 1 - exp( - lambda * pi * r^2)}
  where \eqn{\lambda}{lambda} is the intensity
  (expected number of points per unit area).
  Deviations between the empirical and theoretical \eqn{G} curves
  may suggest spatial clustering or spatial regularity.

  This algorithm estimates the nearest neighbour distance distribution
  function \eqn{G}
  from the point pattern \code{X}. It assumes that \code{X} can be treated
  as a realisation of a stationary (spatially homogeneous) 
  random spatial point process in the plane, observed through
  a bounded window.
  The window (which is specified in \code{X} as \code{X$window})
  may have arbitrary shape. 

  The argument \code{X} is interpreted as a point pattern object 
  (of class \code{"ppp"}, see \code{\link{ppp.object}}) and can
  be supplied in any of the formats recognised
  by \code{\link{as.ppp}()}. 

  The estimation of \eqn{G} is hampered by edge effects arising from 
  the unobservability of points of the random pattern outside the window. 
  An edge correction is needed to reduce bias (Baddeley, 1998; Ripley, 1988). 
  The edge corrections implemented here are the border method or
  ``\emph{reduced sample}'' estimator, the spatial Kaplan-Meier estimator
  (Baddeley and Gill, 1997) and the Hanisch estimator (Hanisch, 1984).

  The argument \code{r} is the vector of values for the
  distance \eqn{r} at which \eqn{G(r)} should be evaluated. 
  It is also used to determine the breakpoints
  (in the sense of \code{\link{hist}})
  for the computation of histograms of distances. The 
  estimators are computed from histogram counts. 
  This introduces a discretisation
  error which is controlled by the fineness of the breakpoints.

  First-time users would be strongly advised not to specify \code{r}.
  However, if it is specified, \code{r} must satisfy \code{r[1] = 0}, 
  and \code{max(r)} must be larger than the radius of the largest disc 
  contained in the window. Furthermore, the successive entries of \code{r}
  must be finely spaced.

  The algorithm also returns an estimate of the hazard rate function, 
  \eqn{\lambda(r)}{lambda(r)}, of \eqn{G(r)}. The hazard rate is
  defined as the derivative
  \deqn{\lambda(r) = - \frac{d}{dr} \log (1 - G(r))}{%
    lambda(r) = - (d/dr) log(1 - G(r))}
  This estimate should be used with caution as \eqn{G} is not necessarily
  differentiable.

  The naive empirical distribution of distances from each point of
  the pattern \code{X} to the nearest other point of the pattern, 
  is a biased estimate of \eqn{G}. However it is sometimes useful.
  It can be returned by the algorithm, by selecting \code{correction="none"}.
  Care should be taken not to use the uncorrected
  empirical \eqn{G} as if it were an unbiased estimator of  \eqn{G}.

  To simply compute the nearest neighbour distance for each point in the
  pattern, use \code{\link{nndist}}. To determine which point is the
  nearest neighbour of a given point, use \code{\link{nnwhich}}.
}
\references{
    Baddeley, A.J. Spatial sampling and censoring.
     In O.E. Barndorff-Nielsen, W.S. Kendall and
     M.N.M. van Lieshout (eds) 
     \emph{Stochastic Geometry: Likelihood and Computation}.
     Chapman and Hall, 1998.
     Chapter 2, pages 37-78.
  
  Baddeley, A.J. and Gill, R.D.
     Kaplan-Meier estimators of interpoint distance
		distributions for spatial point processes.
     \emph{Annals of Statistics} \bold{25} (1997) 263-292.

  Cressie, N.A.C. \emph{Statistics for spatial data}.
    John Wiley and Sons, 1991.

  Diggle, P.J. \emph{Statistical analysis of spatial point patterns}.
    Academic Press, 1983.

  Hanisch, K.-H. (1984) Some remarks on estimators of the distribution
  function of nearest-neighbour distance in stationary spatial point
  patterns. \emph{Mathematische Operationsforschung und Statistik,
    series Statistics} \bold{15}, 409--412.
  
  Ripley, B.D. \emph{Statistical inference for spatial processes}.
  Cambridge University Press, 1988.

  Stoyan, D, Kendall, W.S. and Mecke, J.
  \emph{Stochastic geometry and its applications}.
  2nd edition. Springer Verlag, 1995.
}
\section{Warnings}{
  The function \eqn{G} does not necessarily have a density. 
  Any valid c.d.f. may appear as the nearest neighbour distance
  distribution function of a stationary point process.

  The reduced sample estimator of \eqn{G} is pointwise approximately 
  unbiased, but need not be a valid distribution function; it may 
  not be a nondecreasing function of \eqn{r}. Its range is always 
  within \eqn{[0,1]}.

  The spatial Kaplan-Meier estimator of \eqn{G} is always nondecreasing
  but its maximum value may be less than \eqn{1}.
}
\seealso{
  \code{\link{nndist}},
  \code{\link{nnwhich}},
  \code{\link{Fest}},
  \code{\link{Jest}},
  \code{\link{Kest}},
  \code{\link{km.rs}},
  \code{\link{reduced.sample}},
  \code{\link{kaplan.meier}}
}
\examples{
  data(cells)
  G <- Gest(cells)
  plot(G)

  # P-P style plot
  plot(G, cbind(km,theo) ~ theo)

  # the empirical G is below the Poisson G,
  # indicating an inhibited pattern

  \dontrun{
     plot(G, . ~ r)
     plot(G, . ~ theo)
     plot(G, asin(sqrt(.)) ~ asin(sqrt(theo)))
  }
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{nonparametric}
