\name{nnorient}
\alias{nnorient}
\title{
  Nearest Neighbour Orientation Distribution
}
\description{
  Computes the distribution of the orientation of the vectors
  from each point to its nearest neighbour.
}
\usage{
  nnorient(X, \dots, cumulative = FALSE, correction, k = 1,
                     unit = c("degree", "radian"),
                     domain = NULL, ratio = FALSE)
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"}).
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
    Options are \code{"none"}, \code{"bord.modif"}, 
    \code{"good"} and \code{"best"}.
    Alternatively \code{correction="all"} selects all options.
  }
  \item{k}{
    Integer. The \eqn{k}th nearest neighbour will be used.
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
  This algorithm considers each point in the pattern \code{X}
  and finds its nearest neighbour (or \eqn{k}th nearest neighour).
  The \emph{direction} of the arrow joining the data point to its neighbour
  is measured, as an angle in degrees or radians,
  anticlockwise from the \eqn{x} axis.

  If \code{cumulative=FALSE} (the default),
  a kernel estimate of the probability density of the angles
  is calculated using \code{\link{circdensity}}.
  This is the function \eqn{\vartheta(\phi)}{theta(phi)} defined 
  in Illian et al (2008), equation (4.5.3), page 253.

  If \code{cumulative=TRUE}, then the cumulative distribution
  function of these angles is calculated.

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
  Illian, J., Penttinen, A., Stoyan, H. and Stoyan, D. (2008)
  \emph{Statistical Analysis and Modelling of Spatial Point Patterns.}
  Wiley.
}
\seealso{
  \code{\link{pairorient}}
}
\examples{
  rose(nnorient(redwood, adjust=0.6), col="grey")
  plot(CDF <- nnorient(redwood, cumulative=TRUE))
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}
  
}
\keyword{spatial}
\keyword{nonparametric}
