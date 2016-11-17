\name{qqnorm.im}
\alias{qqnorm.im}
\title{
  Normal Quantile Plot of Pixel Values
}
\description{
  Generates a Normal Quantile-Quantile plot of the pixel values
  of a pixel image.
}
\usage{
  \method{qqnorm}{im}(y, \dots)
}
\arguments{
  \item{y}{
    A pixel image (object of class \code{"im"}).
  }
  \item{\dots}{
    Additional arguments passed to
    \code{\link[stats]{qqnorm.default}}.
  }
}
\details{
  The pixel values of \code{y} are extracted and a normal Q-Q plot
  is generated using \code{\link[stats]{qqnorm.default}}.
}
\value{
   As described for \code{\link[stats]{qqnorm.default}}.
}
\author{
  \adrian
}
\seealso{
  \code{\link{hist.im}}
}
\examples{
  qqnorm(bei.extra$grad, type="l")
}
\keyword{spatial}
\keyword{hplot}
