\name{smooth.msr}
\alias{smooth.msr}
\alias{Smooth.msr}
\title{
  Smooth a Signed or Vector-Valued Measure
}
\description{
  Apply kernel smoothing to a signed measure or vector-valued measure.
}
\usage{
 smooth.msr(X, ...)

 \method{Smooth}{msr}(X, ...)
}
\arguments{
  \item{X}{
    Object of class \code{"msr"} representing a
    signed measure or vector-valued measure.
  }
  \item{\dots}{
    Arguments passed to \code{\link{density.ppp}} controlling the
    smoothing bandwidth and the pixel resolution.
  }
}
\details{
  This function applies kernel smoothing to a signed measure or
  vector-valued measure \code{X}. The Gaussian kernel is used.

  The object \code{X} would typically have been created by
  \code{\link{residuals.ppm}} or \code{\link{msr}}.
}
\value{
  For signed measures, a pixel image (object of class \code{"im"}).
  For vector-valued measures, a list of pixel images; the list also
  belongs to the class \code{"listof"} so that it can be printed and plotted.
}
\references{
  Baddeley, A., Turner, R., Moller, J. and Hazelton, M. (2005)
  Residual analysis for spatial point processes.
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{67}, 617--666.

  Baddeley, A., Moller, J. and Pakes, A.G. (2008) 
  Properties of residuals for spatial point processes.
  \emph{Annals of the Institute of Statistical Mathematics}
  \bold{60}, 627--649.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\seealso{
  \code{\link{Smooth}},
  \code{\link{msr}},
  \code{\link{plot.msr}}
}
\examples{
   example(msr)
   plot(Smooth(rp))
   plot(Smooth(rs))
}
\keyword{spatial}
\keyword{models}
