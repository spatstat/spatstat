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

The result also has an attribute \code{"camefrom"} indicating
the provenance of each line in the resulting network.
For example \code{camefrom[3]=2} means that the third line segment
in the result is a piece of the second segment of \code{X}.
}
\author{
  \spatstatAuthors.
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
  L <- as.linnet(B)

  # check validity
  L
  plot(L)
  text(vertices(L), labels=vertexdegree(L))

  # show the pieces that came from original segment number 1
  S <- as.psp(L)
  (camefrom <- attr(L, "camefrom"))
  parts <- which(camefrom == 1)
  plot(S[parts], add=TRUE, col="green", lwd=2)
}
\keyword{spatial}
\keyword{manip}
