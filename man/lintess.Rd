\name{lintess}
\alias{lintess}
\title{
  Tessellation on a Linear Network
}
\description{
  Create a tessellation on a linear network.
}
\usage{
lintess(L, df, marks=NULL)
}
\arguments{
  \item{L}{
    Linear network (object of class \code{"linnet"}).
  }
  \item{df}{
    Data frame of local coordinates for the pieces that make up the 
    tiles of the tessellation. See Details.
  }
  \item{marks}{
    Vector or data frame of 
    marks associated with the tiles of the tessellation.
  }
}
\details{
  A tessellation on a linear network \code{L} is a partition of the
  network into non-overlapping pieces (tiles). Each tile consists of one
  or more line segments which are subsets of the line segments making up
  the network. A tile can consist of several disjoint pieces.
  
  The data frame \code{df} should have columns named
  \code{seg}, \code{t0}, \code{t1} and \code{tile}.
  Any additional columns will be ignored.

  Each row of the data frame specifies one sub-segment of the network
  and allocates it to a particular tile.
  
  The \code{seg} column specifies which line segment of the network
  contains the sub-segment. Values of \code{seg} are integer indices
  for the segments in \code{as.psp(L)}.

  The \code{t0} and \code{t1} columns specify the start and end points
  of the sub-segment. They should be numeric values between 0 and 1
  inclusive, where the values 0 and 1 representing the network vertices
  that are joined by this network segment.

  The \code{tile} column specifies which tile of the tessellation
  includes this sub-segment. It will be coerced to a factor and its
  levels will be the names of the tiles.

  If \code{df} is missing or \code{NULL}, the result is a tessellation
  with only one tile, consisting of the entire network \code{L}.

  Additional data called \emph{marks} may be associated with
  each tile of the tessellation. The argument \code{marks} should be
  a vector with one entry for each tile (that is, one entry for each
  level of \code{df$tile}) or a data frame with one row for each tile.
  In general \code{df} and \code{marks} will have different numbers of rows.
}
\value{
  An object of class \code{"lintess"}.
  There are methods for \code{print}, \code{plot} and
  \code{summary} for this object.
}
\author{
  \adrian and Greg McSwiggan.
}
\seealso{
  \code{\link{linnet}} for linear networks.

  \code{\link{plot.lintess}} for plotting.

  \code{\link{divide.linnet}} to make a tessellation demarcated by
  given points.

  \code{\link{lineardirichlet}} to create the Dirichlet-Voronoi
  tessellation from a point pattern on a linear network.
  
  \code{\link{as.linfun.lintess}}, \code{\link{as.linnet.lintess}} and
  \code{\link{as.linim}} to convert to other classes.

  \code{\link{tile.lengths}} to compute the length of each tile
  in the tessellation.
  
  The undocumented methods \code{Window.lintess} and
  \code{as.owin.lintess} extract the spatial window.
}
\examples{
   # tessellation consisting of one tile for each existing segment
   ns <- nsegments(simplenet)
   df <- data.frame(seg=1:ns, t0=0, t1=1, tile=letters[1:ns])
   u <- lintess(simplenet, df)
   u
   plot(u)
   S <- as.psp(simplenet)
   marks(u) <- data.frame(len=lengths_psp(S), ang=angles.psp(S))
   u
   plot(u)
}
\keyword{spatial}
\keyword{datagen}

