\name{boundingcircle}
\alias{boundingradius}
\alias{boundingradius.owin}
\alias{boundingradius.ppp}
\alias{boundingcentre}
\alias{boundingcircle}
\alias{boundingcentre.owin}
\alias{boundingcircle.owin}
\alias{boundingcentre.ppp}
\alias{boundingcircle.ppp}
\title{
  Smallest Enclosing Circle
}
\description{
  Find the smallest circle enclosing a spatial window
  or other object. Return its radius, or the location of its centre,
  or the circle itself.
}
\usage{
boundingradius(x, \dots)

boundingcentre(x, \dots)

boundingcircle(x, \dots)

\method{boundingradius}{owin}(x, \dots)

\method{boundingcentre}{owin}(x, \dots)

\method{boundingcircle}{owin}(x, \dots)

\method{boundingradius}{ppp}(x, \dots)

\method{boundingcentre}{ppp}(x, \dots)

\method{boundingcircle}{ppp}(x, \dots)
}
\arguments{
  \item{x}{
    A window (object of class \code{"owin"}), or another spatial object. 
  }
  \item{\dots}{
    Arguments passed to \code{\link{as.mask}} to determine the
    pixel resolution for the calculation.
  }
}
\details{
  The \code{boundingcircle} of a spatial region \eqn{W} is the smallest circle
  that contains \eqn{W}. The \code{boundingradius} is the radius of this circle,
  and the \code{boundingcentre} is the centre of the circle.
  
  The functions \code{boundingcircle}, \code{boundingcentre} and
  \code{boundingradius} are generic. There are methods
  for objects of class \code{"owin"}, \code{"ppp"} and \code{"linnet"}.
}
\value{
  The result of \code{boundingradius} is a single numeric value.

  The result of \code{boundingcentre} is a point pattern containing a
  single point.

  The result of \code{boundingcircle} is a window representing the
  boundingcircle.
}
\author{
  \adrian 
}
\seealso{
 \code{\link{boundingradius.linnet}}
}
\examples{
  boundingradius(letterR)

  plot(grow.rectangle(Frame(letterR), 0.2), main="", type="n")
  plot(letterR, add=TRUE, col="grey")
  plot(boundingcircle(letterR), add=TRUE, border="green", lwd=2)
  plot(boundingcentre(letterR), pch="+", cex=2, col="blue", add=TRUE)

  X <- runifpoint(5)
  plot(X)
  plot(boundingcircle(X), add=TRUE)
  plot(boundingcentre(X), pch="+", cex=2, col="blue", add=TRUE)
}
\keyword{spatial}
\keyword{math}
