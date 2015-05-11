\name{adaptive.density}
\alias{adaptive.density}
\title{Intensity Estimate of Point Pattern Using Tessellation}
\description{
  Computes an adaptive estimate of the intensity function of a point
  pattern.
}
\usage{
adaptive.density(X, f = 0.1, ..., nrep = 1, verbose=TRUE)
}
\arguments{
  \item{X}{Point pattern dataset (object of class \code{"ppp"}).}
  \item{f}{Fraction (between 0 and 1 inclusive) of the data points that will be
    removed from the data and used to determine a tessellation for the
    intensity estimate.
  }
  \item{\dots}{Arguments passed to \code{\link{as.im}} determining the
    pixel resolution of the result.
  }
  \item{nrep}{Number of independent repetitions of the randomised
    procedure.}
  \item{verbose}{Logical value indicating whether to print
    progress reports.}
}
\details{
  This function is an alternative to \code{\link{density.ppp}}. It
  computes an estimate of the intensity function of a point pattern
  dataset. The result is a pixel image giving the estimated intensity,

  If \code{f=1}, the Voronoi estimate (Barr and Schoenberg, 2010)
  is computed: the point pattern \code{X} is used to construct
  a Voronoi/Dirichlet tessellation (see \code{\link{dirichlet}});
  the areas of the Dirichlet tiles are computed; the estimated intensity
  in each tile is the reciprocal of the tile area.

  If \code{f=0}, the intensity estimate at every location is
  equal to the average intensity (number of points divided by window area).
  
  If \code{f} is strictly between 0 and 1, 
  the dataset \code{X} is randomly split into two patterns \code{A} and
  \code{B} containing a fraction \code{f} and \code{1-f}, respectively,
  of the original data. The subpattern \code{A} is used to construct a
  Dirichlet tessellation, while the subpattern
  \code{B} is retained for counting. For each tile of the Dirichlet
  tessellation, we count the number of points of \code{B} falling in the
  tile, and divide by the area of the same tile, to obtain an estimate
  of the intensity of the pattern \code{B} in the tile.
  This estimate is divided by \code{1-f} to obtain an estimate
  of the intensity of \code{X} in the tile. The result is a pixel image
  of intensity estimates which are constant on each tile of the tessellation.

  If \code{nrep} is greater than 1, this randomised procedure is
  repeated \code{nrep} times, and the results are averaged.

  This technique has been used by Ogata et al. (2003), Ogata (2004)
  and Baddeley (2007).
}
\value{
  A pixel image (object of class \code{"im"}) whose values are
  estimates of the intensity of \code{X}.
}
\seealso{
  \code{\link{density.ppp}},
  \code{\link{dirichlet}},
  \code{\link{im.object}}.
}
\references{
  Baddeley, A. (2007)
  Validation of statistical models for spatial point patterns.
  In J.G. Babu and E.D. Feigelson (eds.)
  \emph{SCMA IV: Statistical Challenges in Modern Astronomy IV},
  volume 317 of Astronomical Society of the Pacific Conference Series,
  San Francisco, California USA, 2007. Pages 22--38.

  Barr, C., and Schoenberg, F.P. (2010).
  On the Voronoi estimator for the intensity of an inhomogeneous
  planar Poisson process. \emph{Biometrika} \bold{97} (4), 977--984.

  Ogata, Y. (2004)
  Space-time model for regional seismicity and detection of crustal
  stress changes.
  \emph{Journal of Geophysical Research}, \bold{109}, 2004.

  Ogata, Y., Katsura, K. and Tanemura, M. (2003).
  Modelling heterogeneous space-time occurrences of earthquake and its
  residual analysis.
  \emph{Applied Statistics} \bold{52} 499--509.
}

\examples{
  \dontrun{
  data(nztrees)
  plot(adaptive.density(nztrees))
  }
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
