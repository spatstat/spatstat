\name{rthinclumps}
\alias{rthinclumps}
\title{Random Thinning of Clumps}
\description{
  Finds the topologically-connected clumps of a spatial region
  and randomly deletes some of the clumps.
}
\usage{
rthinclumps(W, p, \dots)
}
\arguments{
  \item{W}{
    Window (object of class \code{"owin"}
    or pixel image (object of class \code{"im"}).
  }
  \item{p}{
    Probability of \emph{retaining} each clump.
    A single number between 0 and 1.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{connected.im}}
    or \code{\link{connected.owin}} to determine the connected clumps.
  }
}
\details{
  The argument \code{W} specifies a region of space, typically
  consisting of several clumps that are not connected to each other.
  The algorithm randomly deletes or retains each clump. The fate of
  each clump is independent of other clumps.

  If \code{W} is a spatial window (class \code{"owin"}) then
  it will be divided into clumps using 
  \code{\link{connected.owin}}. Each clump will either be
  retained (with probability \code{p}) or deleted in its entirety
  (with probability \code{1-p}). 

  If \code{W} is a pixel image (class \code{"im"}) then its domain
  will be divided into clumps using \code{\link{connected.im}}. The default
  behaviour depends on the type of pixel values. If the pixel values are
  logical, then the spatial region will be taken to consist of all
  pixels whose value is \code{TRUE}. Otherwise, the spatial region is
  taken to consist of all pixels whose value is defined (i.e. not equal
  to \code{NA}). This behaviour can be changed using the argument
  \code{background} passed to \code{\link{connected.im}}.

  The result is a window comprising all the clumps that were retained.
}
\value{
  Window (object of class \code{"owin"}).
}
\author{
  \adrian.
}
\seealso{
  \code{\link{rthin}} for thinning other kinds of objects.
}
\examples{
  A <- (distmap(cells) < 0.06)
  opa <- par(mfrow=c(1,2))
  plot(A)
  plot(rthinclumps(A, 0.5))
  par(opa)
}
\keyword{spatial}
\keyword{datagen}
\keyword{manip}

