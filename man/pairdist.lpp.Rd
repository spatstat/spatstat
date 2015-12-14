\name{pairdist.lpp}
\alias{pairdist.lpp}
\title{
  Pairwise shortest-path distances between points on
  a linear network
}
\description{
  Given a pattern of points on a linear network, compute the
  matrix of distances between all pairs of points, measuring
  distance by the shortest path in the network.
}
\usage{
\method{pairdist}{lpp}(X, ..., method="C")
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
  Given a pattern of points on a linear network, this function computes the
  matrix of distances between all pairs of points, measuring
  distance by the shortest path in the network.
  
  If \code{method="C"} the distances are computed using
  code in the C language. If \code{method="interpreted"} then the
  computation is performed using interpreted \R code. The \R code is
  much slower, but is provided for checking purposes.
}
\value{
  A symmetric matrix.
}
\author{
  Ang Qi Wei \email{aqw07398@hotmail.com} and
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
}
\seealso{
  \code{\link{lpp}}
}
\examples{
   example(lpp)
   pairdist(X)
}
\keyword{spatial}
