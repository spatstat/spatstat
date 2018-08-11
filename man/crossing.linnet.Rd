\name{crossing.linnet}
\alias{crossing.linnet}
\title{
  Crossing Points between Linear Network and Other Lines
}
\description{
  Find all the crossing-points between
  a linear network and another pattern of lines or line segments.
}
\usage{
crossing.linnet(X, Y)
}
\arguments{
  \item{X}{
    Linear network (object of class \code{"linnet"}).
  }
  \item{Y}{
    A linear network, or a spatial pattern of line segments (class \code{"psp"})
    or infinite lines (class \code{"infline"}).
  }
}
\details{
  All crossing-points between \code{X} and \code{Y}
  are determined. The result is a point pattern on the network \code{X}.
}
\value{
  Point pattern on a linear network (object of class \code{"lpp"}).
}
\author{
  \adrian.
}
\seealso{
  \code{\link{crossing.psp}}
%  \code{\link{chop.linnet}}
}
\examples{
   plot(simplenet, main="")
   L <- infline(p=runif(3), theta=runif(3, max=pi/2))
   plot(L, col="red")
   Y <- crossing.linnet(simplenet, L)
   plot(Y, add=TRUE, cols="blue")
}
\keyword{spatial}
\keyword{manip}
