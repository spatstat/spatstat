\name{Kmeasure}
\alias{Kmeasure}
\title{Reduced Second Moment Measure}
\description{
Estimates the reduced second moment measure \eqn{\kappa}{Kappa} 
from a point pattern in a window of arbitrary shape.
}
\usage{
  Kmeasure(X, sigma, edge=TRUE, ..., varcov=NULL)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of \eqn{\kappa}{Kappa} will be computed.
    An object of class \code{"ppp"}, or data
    in any format acceptable to \code{\link{as.ppp}()}.
    }
    \item{sigma}{
      Standard deviation \eqn{\sigma}{sigma} of the Gaussian
      smoothing kernel. Incompatible with \code{varcov}.
    }
    \item{edge}{
      logical value indicating whether an edge correction
      should be applied.
    }
    \item{\dots}{Ignored.}
    \item{varcov}{
      Variance-covariance matrix of the Gaussian smoothing kernel.
      Incompatible with \code{sigma}.
    }
}
\value{
  A real-valued pixel image (an object of class \code{"im"},
  see \code{\link{im.object}}) whose pixel values are estimates
  of the density of the reduced second moment measure
  at each location.
}
\details{
  Given a point pattern dataset, 
  this command computes an estimate of the reduced second moment
  measure \eqn{\kappa}{Kappa} of the point process.
  The result is a pixel image whose pixel values are estimates of
  the density of the reduced second moment measure. 

  The reduced second moment measure \eqn{\kappa}{Kappa}
  can be regarded as a generalisation of the more familiar
  \eqn{K}-function.
  An estimate of \eqn{\kappa}{Kappa} derived from a spatial point
  pattern dataset can be useful in exploratory data analysis.
  Its advantage over the \eqn{K}-function is that it is also sensitive
  to anisotropy and directional effects. 
  
  In a nutshell, the command \code{Kmeasure} computes a smoothed version
  of the \emph{Fry plot}. 
  As explained under \code{\link{fryplot}}, the Fry plot is a scatterplot of the
  vectors joining all pairs of points in the pattern.
  The reduced second moment measure is (essentially) defined as
  the average of the Fry plot over different realisations of the point
  process. The command \code{Kmeasure} effectively smooths the Fry plot
  of a dataset to obtain an estimate of the reduced second moment measure.

  In formal terms, the reduced second moment measure \eqn{\kappa}{Kappa}
  of a stationary point process \eqn{X} is a measure defined on the
  two-dimensional plane such that,
  for a `typical' point \eqn{x} of the process, 
  the expected number of other points \eqn{y} of the process
  such that the vector \eqn{y - x} lies in a region \eqn{A},
  equals \eqn{\lambda \kappa(A)}{lambda * Kappa(A)}.
  Here \eqn{\lambda}{lambda}
  is the intensity of the process,
  i.e. the expected number of points of \eqn{X} per unit area.

  The \eqn{K}-function is a special case. The function value \eqn{K(t)} is
  the value of the reduced second moment measure
  for the disc of radius \eqn{t} centred at the origin; that is,
  \eqn{K(t) = \kappa(b(0,t))}{K(t) = Kappa(b(0,t))}.
  
  The command \code{Kmeasure} computes an estimate of \eqn{\kappa}{Kappa}
  from a point pattern dataset \code{X},
  which is assumed to be a realisation of a stationary point process,
  observed inside a known, bounded window. Marks are ignored.

  The algorithm approximates the point pattern and its window by binary pixel
  images, introduces a Gaussian smoothing kernel
  and uses the Fast Fourier Transform \code{\link{fft}}
  to form a density estimate of \eqn{\kappa}{Kappa}. The calculation
  corresponds to the edge correction known as the ``translation
  correction''.

  The Gaussian smoothing kernel may be specified by either of the
  arguments \code{sigma} or \code{varcov}. If \code{sigma} is a single
  number, this specifies an isotropic Gaussian kernel
  with standard deviation \code{sigma} on each coordinate axis.
  If \code{sigma} is a vector of two numbers, this specifies a Gaussian
  kernel with standard deviation \code{sigma[1]} on the \eqn{x} axis,
  standard deviation \code{sigma[2]} on the \eqn{y} axis, and zero
  correlation between the \eqn{x} and \eqn{y} axes. If \code{varcov} is
  given, this specifies the variance-covariance matrix of the
  Gaussian kernel. There do not seem to be any well-established rules
  for selecting the smoothing kernel in this context.
  
  The density estimate of \eqn{\kappa}{Kappa}
  is returned in the form of a real-valued pixel image.
  Pixel values are estimates of the normalised
  second moment density at the centre of the pixel.
  (The uniform Poisson process would have values identically equal to
  \eqn{1}.)
  The image \code{x} and \code{y}
  coordinates are on the same scale as vector displacements in the
  original point pattern window. The point \code{x=0, y=0} corresponds
  to the `typical point'.
  A peak in the image near \code{(0,0)} suggests clustering;
  a dip in the image near \code{(0,0)} suggests inhibition;
  peaks or dips at other positions suggest possible periodicity.
  
  If desired, the value of \eqn{\kappa(A)}{Kappa(A)} for a region
  \eqn{A} can be estimated by computing the integral of the pixel image
  over the domain \eqn{A}, i.e.\ summing the pixel values and
  multiplying by pixel area, using \code{\link{integral.im}}.
  One possible application is to compute anisotropic counterparts of the
  \eqn{K}-function (in which the disc of radius \eqn{t} is replaced
  by another shape). See Examples.
}
\section{Warning}{
  Some writers use the term \emph{reduced second moment measure}
  when they mean the \eqn{K}-function. This has caused
  confusion. 
 
  As originally defined, the
  reduced second moment measure is a measure, obtained by modifying
  the second moment measure, while the \eqn{K}-function is a function
  obtained by evaluating this measure for discs of increasing radius.
  In \pkg{spatstat}, the \eqn{K}-function is computed by
  \code{\link{Kest}} and the reduced second moment measure is computed
  by \code{Kmeasure}.
}
\references{
  Stoyan, D, Kendall, W.S. and Mecke, J. (1995)
  \emph{Stochastic geometry and its applications}.
  2nd edition. Springer Verlag.

  Stoyan, D. and Stoyan, H. (1994)
  Fractals, random shapes and point fields:
  methods of geometrical statistics.
  John Wiley and Sons.
} 
\seealso{
  \code{\link{Kest}},
  \code{\link{fryplot}},
  \code{\link{spatstat.options}},
  \code{\link{integral.im}},
  \code{\link{im.object}}
}
\examples{
 data(cells)
 plot(Kmeasure(cells, 0.05))
 # shows pronounced dip around origin consistent with strong inhibition
 data(redwood)
 plot(Kmeasure(redwood, 0.03), col=grey(seq(1,0,length=32)))
 # shows peaks at several places, reflecting clustering and ?periodicity
 M <- Kmeasure(cells, 0.05)
 # evaluate measure on a sector
 W <- Window(M)
 ang <- as.im(atan2, W)
 rad <- as.im(function(x,y){sqrt(x^2+y^2)}, W)
 sector <- solutionset(ang > 0 & ang < 1 & rad < 0.6)
 integral.im(M[sector, drop=FALSE])
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{nonparametric}
 
 
