\name{pairorient}
\alias{pairorient}
\title{
  Point Pair Orientation Distribution
}
\description{
  Computes the distribution of the orientation of vectors joining
  pairs of points at a particular range of distances.
}
\usage{
pairorient(X, r1, r2, ..., cumulative=FALSE,
           correction, ratio = FALSE,
           unit=c("degree", "radian"), domain=NULL)
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"}).
  }
  \item{r1,r2}{
    Minimum and maximum values of distance
    to be considered.
  }
  \item{\dots}{
    Arguments passed to \code{\link{circdensity}} to control
    the kernel smoothing, if \code{cumulative=FALSE}.
  }
  \item{cumulative}{
    Logical value specifying whether to estimate the probability density
    (\code{cumulative=FALSE}, the default) or the cumulative
    distribution function (\code{cumulative=TRUE}).
  }
  \item{correction}{
    Character vector specifying edge correction or corrections.
    Options are \code{"none"}, \code{"isotropic"}, \code{"translate"},
    \code{"good"} and \code{"best"}.
    Alternatively \code{correction="all"} selects all options.
  }
  \item{ratio}{
    Logical. 
    If \code{TRUE}, the numerator and denominator of
    each edge-corrected estimate will also be saved,
    for use in analysing replicated point patterns.
  }
  \item{unit}{
    Unit in which the angles should be expressed.
    Either \code{"degree"} or \code{"radian"}.
  }
  \item{domain}{
    Optional window. The first point \eqn{x_i}{x[i]} of each pair of points
    will be constrained to lie in \code{domain}.
  }
}
\details{
  This algorithm considers all pairs of points in the pattern
  \code{X} that lie more than \code{r1} and less than \code{r2}
  units apart. The \emph{direction} of the arrow joining the points
  is measured, as an angle in degrees or radians,
  anticlockwise from the \eqn{x} axis.

  If \code{cumulative=FALSE} (the default),
  a kernel estimate of the probability density of the orientations
  is calculated using \code{\link{circdensity}}.

  If \code{cumulative=TRUE}, then the cumulative distribution
  function of these directions is calculated.
  This is the function \eqn{O_{r1,r2}(\phi)}{O[r1,r2](phi)} defined 
  in Stoyan and Stoyan (1994), equation (14.53), page 271.

  In either case the result can be plotted as a rose diagram by 
  \code{\link{rose}}, or as a function plot by \code{\link{plot.fv}}.
  
  The algorithm gives each observed direction a weight,
  determined by an edge correction, to adjust for the fact that some
  interpoint distances are more likely to be observed than others.
  The choice of edge correction or corrections is determined by the argument
  \code{correction}.

  It is also possible to calculate an estimate of the probability
  density from the cumulative distribution function,
  by numerical differentiation. 
  Use \code{\link{deriv.fv}} with the argument \code{Dperiodic=TRUE}.
}
\value{
  A function value table (object of class \code{"fv"})
  containing the estimates of the probability density or the
  cumulative distribution function of angles,
  in degrees (if \code{unit="degree"})
  or radians (if \code{unit="radian"}).
}
\references{
  Stoyan, D. and Stoyan, H. (1994)
  Fractals, random shapes and point fields:
  methods of geometrical statistics.
  John Wiley and Sons.
}
\seealso{
  \code{\link{Kest}}, \code{\link{Ksector}}, \code{\link{nnorient}}
}
\examples{
  rose(pairorient(redwood, 0.05, 0.15, sigma=8), col="grey")
  plot(CDF <- pairorient(redwood, 0.05, 0.15, cumulative=TRUE))
  plot(f <- deriv(CDF, spar=0.6, Dperiodic=TRUE))
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}
  
}
\keyword{spatial}
\keyword{nonparametric}
