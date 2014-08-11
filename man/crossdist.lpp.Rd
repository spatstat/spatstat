\name{crossdist.lpp} %DontDeclareMethods
\alias{crossdist.lpp}
\title{Pairwise distances between two point patterns on a linear network}
\description{
  Computes the distances between pairs of points
  taken from two different point patterns on
  the same linear network.
}
\usage{
  \method{crossdist}{lpp}(X, Y, \dots, method="C")
}
\arguments{
  \item{X,Y}{
    Point patterns on a linear network (objects of class \code{"lpp"}).
    They must lie on the \emph{same} network.
  }
  \item{\dots}{
    Ignored.
  }
  \item{method}{String specifying which method of calculation to use.
    Values are \code{"C"} and \code{"interpreted"}.
  }
}
\value{
  A matrix whose \code{[i,j]} entry is the distance
  from the \code{i}-th point in \code{X}
  to the \code{j}-th point in \code{Y}.
}
\details{
  Given two point patterns on a linear network,
  this function computes the Euclidean distance from each point
  in the first pattern to each point in the second pattern,
  measuring distance by the shortest path in the network.

  This is a method for the generic function \code{\link{crossdist}}
  for point patterns on a linear network
  (objects of class \code{"lpp"}).

  This function expects two point pattern objects \code{X} and \code{Y}
  on the \emph{same} linear network,
  and returns the matrix whose \code{[i,j]} entry is the
  shortest-path distance from \code{X[i]} to \code{Y[j]}.

  The argument \code{method} is not normally used. It is
  retained only for checking the validity of the software.
  If \code{method = "interpreted"} then the distances are
  computed using interpreted R code only. If \code{method="C"}
  (the default) then C code is used. 
  The C code is much faster.
}
\seealso{
  \code{\link{crossdist}},
  \code{\link{crossdist.ppp}},
  \code{\link{pairdist}},
  \code{\link{nndist}}
}
\examples{
   v <- split(chicago)
   X <- v$cartheft
   Y <- v$burglary
   d <- crossdist(X, Y)
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{spatial}
\keyword{math}

