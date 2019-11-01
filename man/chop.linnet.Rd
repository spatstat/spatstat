\name{chop.linnet}
\alias{chop.linnet}
\title{
  Divide a Linear Network into Tiles Using Infinite Lines
}
\description{
  Given a linear network and a set of infinite lines,
  divide the network into tiles demarcated by the lines.
  The result is a tessellation of the network.
}
\usage{
chop.linnet(X, L)
}
\arguments{
  \item{X}{
    Linear network (object of class \code{"linnet"})
    or data acceptable to \code{\link{as.linnet}}.
  }
  \item{L}{
    Infinite line or lines (object of class \code{"infline"}).
  }
}
\details{
  The first line of \code{L} divides \code{X} into two tiles.
  Subsequent lines divide each of these tiles. The result is a
  tessellation of \code{X}. Tiles are not necessarily connected sets.
}
\value{
  Tessellation on a linear network
  (object of class \code{"lintess"}).
}
\author{
  \adrian.
}
\seealso{
  \code{\link{crossing.linnet}}
}
\examples{
   L <- infline(p=runif(3), theta=runif(3, max=pi/2))
   Y <- chop.linnet(simplenet, L)
   plot(Y, main="")
   plot(L, col="red")
}
\keyword{spatial}
\keyword{manip}
