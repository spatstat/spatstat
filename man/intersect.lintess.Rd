\name{intersect.lintess}
\alias{intersect.lintess}
\title{
  Intersection of Tessellations on a Linear Network
}
\description{
  Yields the intersection (common refinement) of two tessellations
  on a linear network.
}
\usage{
intersect.lintess(X, Y)
}
\arguments{
  \item{X,Y}{
    Tessellations (objects of class \code{"lintess"})
    on the same linear network.
  }
}
\value{
  Another tessellation (object of class \code{"lintess"})
  on the same linear network as \code{X} and \code{Y}.
}
\details{
  \code{X} and \code{Y} should be tessellations defined on the same
  linear network.
  
  Each tile in the resulting tessellation is the intersection of
  a tile of \code{X} with a tile of \code{Y}.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{lintess}},
  \code{\link{divide.linnet}},
  \code{\link{chop.linnet}}
}
\examples{
  X <- divide.linnet(runiflpp(4, simplenet))
  Y <- divide.linnet(runiflpp(3, simplenet))
  opa <- par(mfrow=c(1,3))
  plot(X)
  plot(Y)
  plot(intersect.lintess(X,Y))
  par(opa)
}
\keyword{spatial}
\keyword{manip}
