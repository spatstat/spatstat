\name{tile.lengths}
\alias{tile.lengths}
\title{Compute Lengths of Tiles in a Tessellation on a Network}
\description{
  Computes the length of each tile in a tessellation on a linear network.
}
\usage{
tile.lengths(x)
}
\arguments{
  \item{x}{A tessellation on a linear network
    (object of class \code{"lintess"}).}
}
\details{
  A tessellation on a linear network \code{L} is a partition of the
  network into non-overlapping pieces (tiles). Each tile consists of one
  or more line segments which are subsets of the line segments making up
  the network. A tile can consist of several disjoint pieces.

  This command computes the length of each of the tiles 
  that make up the tessellation \code{x}.
  The result is a numeric vector.
}
\value{
  A numeric vector.
}
\seealso{
  \code{\link{lintess}}
}
\examples{
  X <- runiflpp(5, simplenet)
  A <- lineardirichlet(X)
  plot(A)
  tile.lengths(A)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{manip}
