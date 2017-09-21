\name{connected.ppp}  
\Rdversion{1.1}
\alias{connected.ppp}
\title{
  Connected Components of a Point Pattern
}
\description{
  Finds the topologically-connected components of a point pattern,
  when all pairs of points closer than a threshold distance are joined.
}
\usage{
\method{connected}{ppp}(X, R, \dots)
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
    Other arguments, not recognised by these methods.
  }
}
\details{
  This function can be used to identify clumps of points in a point pattern.

  The function \code{connected} is generic. This is the method for
  point patterns (objects of class \code{"ppp"}).

  The point pattern \code{X} is first converted into an abstract graph
  by joining every pair of points that lie closer than \code{R} units
  apart. Then the connected components of this graph are identified.

  Two points in \code{X} belong to the same connected component if they
  can be reached by a series of steps between points of \code{X},
  each step being shorter than \code{R} units in length.

  The result is a vector of labels for the points of \code{X}
  where all the points in a connected component have the same label.
}
\value{
  A point pattern, equivalent to \code{X} except that the points
  have factor-valued marks, with levels corresponding to the
  connected components.
}
\seealso{
  \code{\link{connected.im}},
  \code{\link{im.object}}, 
  \code{\link{tess}}
}
\examples{
   Y <- connected(redwoodfull, 0.1)
   if(interactive()) {
    plot(Y, cols=1:length(levels(marks(Y))),
         main="connected(redwoodfull, 0.1)")
   }
}
\author{
  \adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{math}
