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
\method{nnwhich}{lpp}(X, ..., k=1, method="C")
}
\arguments{
  \item{X}{
    Point pattern on linear network (object of class \code{"lpp"}).
  }
  \item{method}{
    Optional string determining the method of calculation.
    Either \code{"interpreted"} or \code{"C"}.
  }
  \item{k}{
    Integer, or integer vector. The algorithm will find the
    \code{k}th nearest neighbour. 
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

  Alternatively a matrix with one row for each point in \code{X}
  and one column for each entry of \code{k}.
}
\author{
  \adrian
}
\seealso{
  \code{\link{lpp}}
}
\examples{
   X <- runiflpp(10, simplenet)
   nnwhich(X)
   nnwhich(X, k=2)
}
\keyword{spatial}
