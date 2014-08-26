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
pairorient(X, r1, r2, ..., correction, ratio = FALSE)
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
    Ignored.
  }
  \item{correction}{
    Character vector specifying edge correction or corrections.
    Options are \code{"none"}, \code{"isotropic"}, \code{"translate"},
    \code{"good"} and \code{"best"}.
  }
  \item{ratio}{
    Logical. 
    If \code{TRUE}, the numerator and denominator of
    each edge-corrected estimate will also be saved,
    for use in analysing replicated point patterns.
  }
}
\details{
  This function calculates the point pair orientation distribution
  function \eqn{O_{r1,r2}(\phi)}{O[r1,r2](phi)} defined 
  in Stoyan and Stoyan (1994), equation (14.53), page 271.

  The function considers all pairs of points in the pattern
  \code{X} that lie more than \code{r1} and less than \code{r2}
  units apart. The \emph{direction} of the arrow joining the points
  is measured, as an angle in degrees, anticlockwise from the \eqn{x} axis.
  The result is the cumulative distribution function of these
  directions. 

  In calculating the cumulative distribution function,
  the algorithm gives each observed direction a weight,
  determined by an edge correction, to adjust for the fact that some
  interpoint distances are more likely to be observed than others.
  The choice of edge correction or corrections is determined by the argument
  \code{correction}.

  To calculate the probability density of directions,
  use \code{\link{deriv.fv}} with the argument \code{Dperiodic=TRUE}.
}
\value{
  A function value table (object of class \code{"fv"})
  containing the estimates of the cumulative distribution
  function of angles, in degrees.
}
\references{
  Stoyan, D. and Stoyan, H. (1994)
  Fractals, random shapes and point fields:
  methods of geometrical statistics.
  John Wiley and Sons.
}
\seealso{
  \code{\link{Kest}}, \code{\link{Ksector}}
}
\examples{
  plot(f <- pairorient(redwood, 0.05, 0.15))
  plot(Df <- deriv(f, spar=0.6, Dperiodic=TRUE))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{nonparametric}
