\name{densityVoronoi}
\alias{densityVoronoi}
\alias{densityVoronoi.ppp}
\title{Intensity Estimate of Point Pattern Using Voronoi-Dirichlet Tessellation}
\description{
  Computes an adaptive estimate of the intensity function of a point
  pattern using the Dirichlet-Voronoi tessellation.
}
\usage{
densityVoronoi(X, \dots)

\method{densityVoronoi}{ppp}(X, f = 1, \dots,
                          counting=FALSE,
                          fixed=FALSE,
                          nrep = 1, verbose=TRUE)
}
\arguments{
  \item{X}{Point pattern dataset (object of class \code{"ppp"}).}
  \item{f}{
    Fraction (between 0 and 1 inclusive) of the data points that will be
    used to build a tessellation for the intensity estimate.
  }
  \item{\dots}{Arguments passed to \code{\link{as.im}} determining the
    pixel resolution of the result.
  }
  \item{counting}{
    Logical value specifying the choice of estimation method.
    See Details.
  }
  \item{fixed}{
    Logical. If \code{FALSE} (the default), the data points are independently
    randomly thinned, so the number of data points that are retained
    is random. If \code{TRUE}, the number of data points retained
    is fixed. See Details.
  }
  \item{nrep}{Number of independent repetitions of the randomised
    procedure.}
  \item{verbose}{Logical value indicating whether to print
    progress reports.}
}
\details{
  This function is an alternative to \code{\link{density.ppp}}. It
  computes an estimate of the intensity function of a point pattern
  dataset. The result is a pixel image giving the estimated intensity.

  If \code{f=1} (the default), the Voronoi estimate (Barr and Schoenberg, 2010)
  is computed: the point pattern \code{X} is used to construct
  a Voronoi/Dirichlet tessellation (see \code{\link{dirichlet}});
  the areas of the Dirichlet tiles are computed; the estimated intensity
  in each tile is the reciprocal of the tile area.
  The result is a pixel image
  of intensity estimates which are constant on each tile of the tessellation.

  If \code{f=0}, the intensity estimate at every location is
  equal to the average intensity (number of points divided by window area).
  The result is a pixel image
  of intensity estimates which are constant.
  
  If \code{f} is strictly between 0 and 1,
  the estimation method is applied to a random subset of \code{X}.
  This randomised procedure is repeated \code{nrep} times,
  and the results are averaged.
  The subset is selected as follows:
  \itemize{
    \item
    if \code{fixed=FALSE},
    the dataset \code{X} is randomly
    thinned by deleting or retaining each point independently, with
    probability \code{f} of retaining a point.
    \item
    if \code{fixed=TRUE},
    a random sample of fixed size \code{m} is taken from
    the dataset \code{X}, where \code{m} is the largest integer
    less than or equal to \code{f*n} and \code{n} is the number of
    points in \code{X}.
  }
  Then the intensity estimate is calculated as follows:
  \itemize{
    \item if \code{counting = FALSE} (the default), the thinned pattern
    is used to construct a Dirichlet tessellation and form the
    Voronoi estimate (Barr and Schoenberg, 2010) which is then
    adjusted by a factor \code{1/f} or \code{n/m} as appropriate.
    to obtain an estimate
    of the intensity of \code{X} in the tile.
    \item if \code{counting = TRUE},
    the randomly selected subset \code{A}
    is used to construct a Dirichlet tessellation, while the
    complementary subset \code{B} (consisting of points that were not
    selected in the sample) is used for counting
    to calculate a quadrat count estimate of intensity.
    For each tile of the Dirichlet tessellation formed by \code{A},
    we count the number of points of \code{B} falling in the
    tile, and divide by the area of the same tile, to obtain an estimate
    of the intensity of the pattern \code{B} in the tile.
    This estimate is adjusted by \code{1/(1-f)}
    or \code{n/(n-m)} as appropriate 
    to obtain an estimate of the intensity of \code{X} in the tile.
  }

  Ogata et al. (2003) and Ogata (2004) estimated intensity using the
  Dirichlet-Voronoi tessellation in a modelling context.
  Baddeley (2007) proposed intensity estimation by subsampling
  with \code{0 < f < 1}, and used the  technique described above 
  with \code{fixed=TRUE} and \code{counting=TRUE}.
  Barr and Schoenberg (2010) described and analysed the
  Voronoi estimator (corresponding to \code{f=1}).
  Moradi et al (2019) developed the subsampling technique with
  \code{fixed=FALSE} and \code{counting=FALSE} and called it the
  \emph{smoothed Voronoi estimator}.
}
\value{
  A pixel image (object of class \code{"im"}) whose values are
  estimates of the intensity of \code{X}.
}
\seealso{
  \code{\link{densityVoronoi.lpp}},
  \code{\link{adaptive.density}},
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

  Moradi, M., Cronie, 0., Rubak, E., Lachieze-Rey, R.,
  Mateu, J. and Baddeley, A. (2019)
  Resample-smoothing of Voronoi intensity estimators.
  \emph{Statistics and Computing}, in press.

  Ogata, Y. (2004)
  Space-time model for regional seismicity and detection of crustal
  stress changes.
  \emph{Journal of Geophysical Research}, \bold{109}, 2004.

  Ogata, Y., Katsura, K. and Tanemura, M. (2003).
  Modelling heterogeneous space-time occurrences of earthquakes and its
  residual analysis.
  \emph{Applied Statistics} \bold{52} 499--509.
}
\examples{
  plot(densityVoronoi(nztrees, 1, f=1), main="Voronoi estimate")
  nr <- if(interactive()) 100 else 5
  plot(densityVoronoi(nztrees, f=0.5, nrep=nr), main="smoothed Voronoi estimate")
}
\author{
  \spatstatAuthors and Mehdi Moradi.
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
