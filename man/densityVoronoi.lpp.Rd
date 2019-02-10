\name{densityVoronoi.lpp}
\alias{densityVoronoi.lpp}
\title{Intensity Estimate of Point Pattern on Linear Network Using Voronoi-Dirichlet Tessellation}
\description{
  Computes an adaptive estimate of the intensity function of a point
  pattern on a linear network, using the Dirichlet-Voronoi tessellation
  on the network.
}
\usage{
\method{densityVoronoi}{lpp}(X, f = 1, \dots, nrep = 1, verbose = TRUE)
}
\arguments{
  \item{X}{
    Point pattern on a linear network (object of class \code{"lpp"}).
  }
  \item{f}{
    Fraction (between 0 and 1 inclusive) of the data points that will be
    used to build a tessellation for the intensity estimate.
  }
  \item{\dots}{
    Arguments passed to \code{\link{linim}} determining the
    pixel resolution of the result.
  }
  \item{nrep}{
    Number of independent repetitions of the randomised
    procedure.
  }
  \item{verbose}{
    Logical value indicating whether to print
    progress reports.
  }
}
\details{
  This function is an alternative to \code{\link{density.lpp}}. It
  computes an estimate of the intensity function of a point pattern
  dataset on a linear network.
  The result is a pixel image on the network, giving the estimated intensity.

  This function is a method for the generic \code{\link{densityVoronoi}}
  for the class \code{"lpp"} of point patterns on a linear network.
  
  If \code{f=1} (the default), the Voronoi estimate (Barr and Schoenberg, 2010)
  is computed: the point pattern \code{X} is used to construct
  a Voronoi/Dirichlet tessellation on the network
  (see \code{\link{lineardirichlet}});
  the lengths of the Dirichlet tiles are computed; the estimated intensity
  in each tile is the reciprocal of the tile length.
  The result is a pixel image
  of intensity estimates which are constant on each tile of the tessellation.

  If \code{f=0}, the intensity estimate at every location is
  equal to the average intensity (number of points divided by
  network length). The result is a pixel image
  of intensity estimates which are constant.
  
  If \code{f} is strictly between 0 and 1,
  the smoothed Voronoi estimate (Moradi et al, 2019) is computed.
  The dataset \code{X} is randomly
  thinned by deleting or retaining each point independently, with
  probability \code{f} of retaining a point.
  The thinned pattern
  is used to construct a Dirichlet tessellation and form the
  Voronoi estimate, which is then
  adjusted by a factor \code{1/f}.
  This procedure is repeated \code{nrep} times and the results are
  averaged to obtain the smoothed Voronoi estimate.

  The value \code{f} can be chosen automatically by bandwidth
  selection using \code{\link{bw.voronoi}}.
}
\value{
  Pixel image on a linear network (object of class \code{"linim"}).
}
\references{
  Moradi, M., Cronie, 0., Rubak, E., Lachieze-Rey, R.,
  Mateu, J. and Baddeley, A. (2019)
  Resample-smoothing of Voronoi intensity estimators.
  \emph{Statistics and Computing}, in press.
}
\author{
  \spatstatAuthors and Mehdi Moradi.
}
\seealso{
  \code{\link{densityVoronoi}} is the generic, with a method for
  class \code{"ppp"}.

  \code{\link{lineardirichlet}} computes the Dirichlet-Voronoi
  tessellation on a network.
  
  \code{\link{bw.voronoi}} performs bandwidth selection of the fraction \code{f}.

  See also \code{\link{density.lpp}}.
}
\examples{
   nr <- if(interactive()) 100 else 3
   plot(densityVoronoi(spiders, 0.1, nrep=nr))
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
