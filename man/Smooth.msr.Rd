\name{Smooth.msr}
\alias{Smooth.msr}
\title{
  Smooth a Signed or Vector-Valued Measure
}
\description{
  Apply kernel smoothing to a signed measure or vector-valued measure.
}
\usage{
 \method{Smooth}{msr}(X, ..., drop=TRUE)
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
  \item{drop}{
    Logical. If \code{TRUE} (the default), the result of smoothing
    a scalar-valued measure is a pixel image. If \code{FALSE}, the
    result of smoothing a scalar-valued measure is a list
    containing one pixel image.
  }
}
\details{
  This function applies kernel smoothing to a signed measure or
  vector-valued measure \code{X}. The Gaussian kernel is used.

  The object \code{X} would typically have been created by
  \code{\link{residuals.ppm}} or \code{\link{msr}}.
}
\value{
  A pixel image or a list of pixel images.
  For scalar-valued measures, a pixel image (object of class
  \code{"im"}) provided \code{drop=TRUE}.
  For vector-valued measures (or if \code{drop=FALSE}),
  a list of pixel images; the list also
  belongs to the class \code{"solist"} so that it can be printed and plotted.
}
\references{
  Baddeley, A., Turner, R.,
  \ifelse{latex}{\out{M\o ller}}{Moller}, J. and Hazelton, M. (2005)
  Residual analysis for spatial point processes.
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{67}, 617--666.

  Baddeley, A., \ifelse{latex}{\out{M\o ller}}{Moller}, J. and Pakes, A.G. (2008) 
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
