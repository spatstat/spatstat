\name{divide.linnet}
\alias{divide.linnet}
\title{
  Divide Linear Network at Cut Points
}
\description{
  Make a tessellation of a linear network by dividing it into
  pieces demarcated by the points of a point pattern.
}
\usage{
 divide.linnet(X)
}
\arguments{
  \item{X}{
    Point pattern on a linear network (object of class \code{"lpp"}).
  }
}
\details{
  The points \code{X} are interpreted as dividing the linear network
  \code{L=as.linnet(X)} into separate pieces.

  Two locations on \code{L}
  belong to the same piece if and only if they can be joined by a path
  in \code{L} that does not cross any of the points of \code{X}.

  The result is a
  tessellation of the network (object of class \code{"lintess"})
  representing the division of \code{L} into pieces. 
}
\value{
  A tessellation on a linear network (object of class \code{"lintess"}).
}
\author{
  \spatstatAuthors
  and Greg McSwiggan.
}
\seealso{
  \code{\link{linnet}}, \code{\link{lintess}}.
}
\examples{
  X <- runiflpp(5, simplenet)
  plot(divide.linnet(X))
  plot(X, add=TRUE, pch=16, show.network=FALSE)
}
\keyword{spatial}
\keyword{manip}
