\name{as.data.frame.lintess}
\alias{as.data.frame.lintess}
\title{Convert Network Tessellation to Data Frame}
\description{
  Converts a tessellation on a linear network into a data frame.
}
\usage{
\method{as.data.frame}{lintess}(x, \dots)
}
\arguments{
  \item{x}{
    Tessellation on a linear network (object of class \code{"lintess"}).
  }
  \item{\dots}{Further arguments passed to
    \code{\link[base:as.data.frame]{as.data.frame.default}}
    to determine the row names and other features.
  }
}
\details{
  A tessellation on a linear network is a partition of the
  network into non-overlapping pieces (tiles). Each tile consists of one
  or more line segments which are subsets of the line segments making up
  the network. A tile can consist of several disjoint pieces.

  This function converts the tessellation \code{x} to a data frame.
  Each row of the data frame specifies one sub-segment of the network,
  and allocates it to a particular tile. The data frame has the
  following columns:
  \itemize{
    \item The \code{seg} column specifies which line segment of the network
    contains the sub-segment. Values of \code{seg} are integer indices
    for the network segments in \code{as.psp(as.linnet(x))}.
    \item The \code{t0} and \code{t1} columns specify the start and end points
    of the sub-segment. They are numeric values between 0 and 1
    inclusive, where the values 0 and 1 representing the network vertices
    that are joined by this network segment.
    \item The \code{tile} column specifies which tile of the tessellation
    includes this sub-segment. It is a factor whose 
    levels are the names of the tiles.
  }

  The tessellation may have marks, which are attached to the \emph{tiles}
  of the tessellation. If marks are present, the resulting data frame
  includes columns containing, for each sub-segment,
  the mark value of the corresponding tile.
}
\value{
  A data frame with columns named \code{seg}, \code{t0}, \code{t1},
  \code{tile}, and possibly other columns.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{lintess}}
}
\examples{
  X <- lineardirichlet(runiflpp(3, simplenet))
  marks(X) <- letters[1:3]
  as.data.frame(X)
}
\keyword{spatial}
\keyword{methods}

