\name{nndist.psp}
\alias{nndist.psp}
\title{Nearest neighbour distances between line segments}
\description{
  Computes the distance from each line segment 
  to its nearest neighbour in a line segment pattern.
  Alternatively finds the distance to the second nearest,
  third nearest etc.
}
\usage{
  \method{nndist}{psp}(X, \dots, k=1, method="C")
}
\arguments{
  \item{X}{
    A line segment pattern (object of class \code{"psp"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{k}{
    Integer, or integer vector. The algorithm will compute the distance to the
    \code{k}th nearest neighbour. 
  }
  \item{method}{
    String specifying which method of calculation to use.
    Values are \code{"C"} and \code{"interpreted"}.
    Usually not specified.
  }
}
\value{
  Numeric vector or matrix containing the 
  nearest neighbour distances for each line segment.

  If \code{k = 1} (the default), the return value is a
  numeric vector \code{v} such that \code{v[i]} is the
  nearest neighbour distance for the \code{i}th segment.
  
  If \code{k} is a single integer, then the return value is a
  numeric vector \code{v} such that \code{v[i]} is the
  \code{k}th nearest neighbour distance for the
  \code{i}th segment.

  If \code{k} is a vector, then the return value is a
  matrix \code{m} such that \code{m[i,j]} is the
  \code{k[j]}th nearest neighbour distance for the
  \code{i}th segment.
}
\details{
  This is a method for the generic function \code{\link{nndist}}
  for the class \code{"psp"}.

  If \code{k=1}, this function computes the distance from each line segment 
  to the nearest other line segment in \code{X}.
  In general it computes the distance from each line segment to the
  \code{k}th nearest other line segment.
  The argument \code{k} can also be a vector, and this computation will
  be performed for each value of \code{k}.

  Distances are calculated using the Hausdorff metric. The Hausdorff
  distance between two line segments is the maximum distance
  from any point on one of the segments to the nearest point on
  the other segment.   

  If there are fewer than \code{max(k)+1} line segments in the pattern,
  some of the nearest neighbour distances will be infinite (\code{Inf}).
  
  The argument \code{method} is not normally used. It is
  retained only for checking the validity of the software.
  If \code{method = "interpreted"} then the distances are
  computed using interpreted \R code only. If \code{method="C"}
  (the default) then compiled \code{C} code is used.
  The \code{C} code is somewhat faster.
}
\seealso{
  \code{\link{nndist}},
  \code{\link{nndist.ppp}}
}
\examples{
   L <- psp(runif(10), runif(10), runif(10), runif(10), owin())
   D <- nndist(L)
   D <- nndist(L, k=1:3)
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}

