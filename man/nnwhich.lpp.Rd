\name{nnwhich.lpp}
\alias{nnwhich.lpp}
\title{
  Identify Nearest Neighbours on a Linear Network
}
\description{
  Given a pattern of points on a linear network, identify the
  nearest neighbour for each point, measured
  by the shortest path in the network.
}
\usage{
\method{nnwhich}{lpp}(X, ..., method="C")
}
\arguments{
  \item{X}{
    Point pattern on linear network (object of class \code{"lpp"}).
  }
  \item{method}{
    Optional string determining the method of calculation.
    Either \code{"interpreted"} or \code{"C"}.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  Given a pattern of points on a linear network,
  this function finds the nearest neighbour of each point
  (i.e. for each point it identifies the nearest other point)
  measuring distance by the shortest path in the network.
  
  If \code{method="C"} the task is performed using
  code in the C language. If \code{method="interpreted"} then the
  computation is performed using interpreted \R code. The \R code is
  much slower, but is provided for checking purposes.
}
\value{
  An integer vector, of length equal to the number of points in
  \code{X}, identifying the nearest neighbour of each point.
  If \code{nnwhich(X)[2] = 4} then the nearest neighbour of
  point 2 is point 4.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\seealso{
  \code{\link{lpp}}
}
\examples{
   example(lpp)
   nnwhich(X)
}
\keyword{spatial}
