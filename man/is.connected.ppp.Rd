\name{is.connected.ppp}  
\Rdversion{1.1}
\alias{is.connected.ppp}
\title{
  Determine Whether a Point Pattern is Connected
}
\description{
  Determine whether a point pattern is topologically connected
  when all pairs of points closer than a threshold distance are joined.
}
\usage{
\method{is.connected}{ppp}(X, R, \dots)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
  \item{R}{
    Threshold distance. Pairs of points closer than \code{R} units apart
    will be joined together.
  }
  \item{\dots}{
     Ignored.
  }
}
\details{
  The function \code{is.connected} is generic. This is the method for
  point patterns (objects of class \code{"ppp"}).

  The point pattern \code{X} is first converted into an abstract graph
  by joining every pair of points that lie closer than \code{R} units
  apart. Then the algorithm determines whether this graph is connected.

  That is, the result of \code{is.connected(X)} is \code{TRUE}
  if any point in \code{X} can be reached from any other point,
  by a series of steps between points of \code{X},
  each step being shorter than \code{R} units in length.
}
\value{
  A logical value.
}
\seealso{
  \code{\link{is.connected}}, 
  \code{\link{connected.ppp}}.
}
\examples{
  is.connected(redwoodfull, 0.1)
  is.connected(redwoodfull, 0.2)
}
\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{math}
