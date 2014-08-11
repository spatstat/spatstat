\name{plot.scan.test}
\alias{plot.scan.test}
\alias{as.im.scan.test}
\title{
  Plot Result of Scan Test
}
\description{
  Computes or plots an image showing the
  likelihood ratio test statistic for the scan test,
  or the optimal circle radius. 
}
\usage{
 \method{plot}{scan.test}(x, \dots, what=c("statistic", "radius"),
               do.window = TRUE)

 \method{as.im}{scan.test}(X, \dots, what=c("statistic", "radius"))
}
\arguments{
  \item{x,X}{
    Result of a scan test. An object of class \code{"scan.test"}
    produced by \code{\link{scan.test}}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{plot.im}} to control the appearance
    of the plot.
  }
  \item{what}{
    Character string indicating whether to produce an image of the
    (profile) likelihood ratio test statistic (\code{what="statistic"},
    the default) or an image of the optimal value of circle radius
    (\code{what="radius"}).
  }
  \item{do.window}{
    Logical value indicating whether to plot the original window
    of the data as well.
  }
}
\details{
  These functions extract, and plot, the spatially-varying value
  of the likelihood ratio test statistic which forms the basis of
  the scan test.

  If the test result \code{X} was based on circles of
  the same radius \code{r}, then \code{as.im(X)} is a pixel image
  of the likelihood ratio test statistic as a function of the
  position of the centre of the circle.

  If the test result \code{X} was based on circles of
  several different radii \code{r}, then \code{as.im(X)} is a pixel image
  of the profile (maximum value over all radii \code{r})
  likelihood ratio test statistic as a function of the
  position of the centre of the circle, and
  \code{as.im(X, what="radius")} is a pixel image giving
  for each location \eqn{u} the value of \code{r} which maximised
  the likelihood ratio test statistic at that location.

  The \code{plot} method plots the corresponding image.
}
\value{
  The value of \code{as.im.scan.test} is a pixel image (object of
  class \code{"im"}). The value of \code{plot.scan.test} is \code{NULL}.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{scan.test}},
  \code{\link{scanLRTS}}
}
\examples{
   if(interactive()) {
     a <- scan.test(redwood, seq(0.04, 0.1, by=0.01),
                    method="poisson", nsim=19)
   } else {
     a <- scan.test(redwood, c(0.05, 0.1), method="poisson", nsim=2)
   }
   plot(a)
   as.im(a)
   plot(a, what="radius")
}
\keyword{htest}
\keyword{spatial}
