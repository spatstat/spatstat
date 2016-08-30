\name{as.linnet.psp}
\alias{as.linnet.psp}
\title{
  Convert Line Segment Pattern to Linear Network
}
\description{
  Converts a line segment pattern to a linear network.
}
\usage{
 \method{as.linnet}{psp}(X, \dots, eps, sparse=FALSE)
}
\arguments{
  \item{X}{
    Line segment pattern (object of class \code{"psp"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{eps}{
    Optional. Distance threshold. If two segment endpoints
    are closer than \code{eps} units apart, they will be treated
    as the same point, and will become a single vertex in the
    linear network.
  }
  \item{sparse}{
    Logical value indicating whether to use a sparse matrix
    representation, as explained in \code{\link{linnet}}.
  }
}
\details{
  This command converts any collection of line segments into a linear
  network by guessing the connectivity of the network,
  using the distance threshold \code{eps}.

  If any segments in \code{X} cross over each other, they are first
  cut into pieces using \code{\link{selfcut.psp}}.

  Then any pair of segment endpoints lying closer than \code{eps}
  units apart, is treated as a single vertex. The linear network
  is then constructed using \code{\link{linnet}}.

  It would be wise to check the result by plotting the degree
  of each vertex, as shown in the Examples.

  If \code{X} has marks, then these are stored in
  the resulting linear network \code{Y <- as.linnet(X)},
  and can be extracted as \code{marks(as.psp(Y))} or \code{marks(Y$lines)}.
}
\value{
A linear network (object of class \code{"linnet"}).
}
\author{
\adrian


\rolf

and \ege

}
\seealso{
  \code{\link{linnet}}, \code{\link{selfcut.psp}},
  \code{\link{methods.linnet}}.
}
\examples{
  # make some data
  A <- psp(0.09, 0.55, 0.79, 0.80, window=owin())
  B <- superimpose(A, as.psp(simplenet))

  # convert to a linear network
  D <- as.linnet(B)

  # check validity
  D
  plot(D)
  text(vertices(D), labels=vertexdegree(D))
}
\keyword{spatial}
\keyword{manip}
