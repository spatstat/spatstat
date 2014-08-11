\name{plot.quadrattest}
\alias{plot.quadrattest}
\title{
  Display the result of a quadrat counting test.
}
\description{
  Given the result of a quadrat counting test,
  graphically display the quadrats that were used, the 
  observed and expected counts, and the residual in each quadrat.
}
\usage{
 \method{plot}{quadrattest}(x, ..., textargs=list())
}
\arguments{
  \item{x}{
    Object of class \code{"quadrattest"} containing the result
    of \code{\link{quadrat.test}}.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{plot.tess}} to
    control the display of the quadrats.
  }
  \item{textargs}{
    List of additional arguments passed to
    \code{\link[graphics]{text.default}}
    to control the appearance of the text.
  }
}
\details{
  This is the plot method for objects
  of class \code{"quadrattest"}. Such an object is produced by
  \code{\link{quadrat.test}} and represents the result of
  a \eqn{\chi^2}{chi^2} test for a spatial point pattern.

  The quadrats are first plotted using \code{\link{plot.tess}}.
  Then in each quadrat, the observed and expected counts
  and the Pearson residual are displayed as text using
  \code{\link[graphics]{text.default}}.
  Observed count is displayed at top left; expected count at top right;
  and Pearson residual at bottom. 
}
\value{
  Null.
}
\examples{
   plot(quadrat.test(swedishpines, 3))
}
\seealso{
  \code{\link{quadrat.test}},
  \code{\link{plot.tess}},
  \code{\link[graphics]{text.default}},
  \code{\link{plot.quadratcount}}
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{htest}
\keyword{hplot}
