\name{nnwhich}
\alias{nnwhich}
\alias{nnwhich.ppp}
\alias{nnwhich.default}
\title{Nearest neighbour}
\description{
  Finds the nearest neighbour of each point in a point pattern.
}
\usage{
  nnwhich(X, \dots)
  \method{nnwhich}{ppp}(X, \dots, k=1, by=NULL, method="C")
  \method{nnwhich}{default}(X, Y=NULL, \dots, k=1, by=NULL, method="C")
}
\arguments{
  \item{X,Y}{
    Arguments specifying the locations of
    a set of points.
    For \code{nnwhich.ppp}, the argument \code{X} should be a point
    pattern (object of class \code{"ppp"}).
    For \code{nnwhich.default}, typically \code{X} and \code{Y} would be
    numeric vectors of equal length. Alternatively \code{Y} may be
    omitted and \code{X} may be
    a list with two components \code{x} and \code{y},
    or a matrix with two columns.
  }
  \item{\dots}{
    Ignored by \code{nnwhich.ppp}
    and \code{nnwhich.default}.
  }
  \item{k}{
    Integer, or integer vector. The algorithm will compute the distance to the
    \code{k}th nearest neighbour. 
  }
  \item{by}{
    Optional. A factor, which separates \code{X} into groups.
    The algorithm will find the nearest neighbour in each group. 
  }
  \item{method}{String specifying which method of calculation to use.
    Values are \code{"C"} and \code{"interpreted"}.
  }
}
\value{
  Numeric vector or matrix giving, for each point,
  the index of its nearest neighbour (or \code{k}th nearest neighbour).

  If \code{k = 1} (the default), the return value is a
  numeric vector \code{v} giving the indices of the nearest neighbours
  (the nearest neighbout of the \code{i}th point is
  the \code{j}th point where \code{j = v[i]}).
  
  If \code{k} is a single integer, then the return value is a
  numeric vector giving the indices of the
  \code{k}th nearest neighbours.

  If \code{k} is a vector, then the return value is a
  matrix \code{m} such that \code{m[i,j]} is the
  index of the \code{k[j]}th nearest neighbour for the
  \code{i}th data point.

  If the argument \code{by} is given, then the result is a data frame
  containing the indices described above, from each point of \code{X},
  to the nearest point in each subset of \code{X}
  defined by the factor \code{by}.
}
\details{
  For each point in the given point pattern, this function finds
  its nearest neighbour (the nearest other point of the pattern).
  By default it returns a vector giving, for each point,
  the index of the point's
  nearest neighbour. If \code{k} is specified, the algorithm finds
  each point's \code{k}th nearest neighbour.

  The function \code{nnwhich} is generic, with
  method for point patterns (objects of class \code{"ppp"})
  and a default method which are described here, as well as a method for
  three-dimensional point patterns (objects of class \code{"pp3"},
  described in \code{\link{nnwhich.pp3}}.

  The method \code{nnwhich.ppp} expects a single
  point pattern argument \code{X}.
  The default method expects that \code{X} and \code{Y} will determine
  the coordinates of a set of points. Typically \code{X} and
  \code{Y} would be numeric vectors of equal length. Alternatively
  \code{Y} may be omitted and \code{X} may be a list with two components
  named \code{x} and \code{y}, or a matrix or data frame with two columns.
  
  The argument \code{k} may be a single integer, or an integer vector.
  If it is a vector, then the \eqn{k}th nearest neighbour distances are
  computed for each value of \eqn{k} specified in the vector.

  If the argument \code{by} is given, it should be a \code{factor},
  of length equal to the number of points in \code{X}.
  This factor effectively partitions \code{X} into subsets,
  each subset associated with one of the levels of \code{X}.
  The algorithm will then find, for each point of \code{X},
  the nearest neighbour \emph{in each subset}.

  If there are no points (if \code{x} has length zero)
  a numeric vector of length zero is returned.
  If there is only one point (if \code{x} has length 1),
  then the nearest neighbour is undefined, and a value of \code{NA}
  is returned. In general if the number of points is less than or equal
  to \code{k}, then a vector of \code{NA}'s is returned.

  The argument \code{method} is not normally used. It is
  retained only for checking the validity of the software.
  If \code{method = "interpreted"} then the distances are
  computed using interpreted R code only. If \code{method="C"}
  (the default) then C code is used. 
  The C code is faster by two to three orders of magnitude
  and uses much less memory.
  
  To evaluate the \emph{distance} between a point and its nearest
  neighbour, use \code{\link{nndist}}.

  To find the nearest neighbours from one point pattern
  to another point pattern, use \code{\link{nncross}}.
}
\section{Nearest neighbours of each type}{
  If \code{X} is a multitype point pattern 
  and \code{by=marks(X)}, then the algorithm will find,
  for each point of \code{X}, the nearest neighbour
  of each type. See the Examples. 
}
\section{Warnings}{
  A value of \code{NA} is returned if there is only one point
  in the point pattern. 
}
\seealso{
  \code{\link{nndist}},
  \code{\link{nncross}}
}
\examples{
   data(cells)
   plot(cells)
   m <- nnwhich(cells)
   m2 <- nnwhich(cells, k=2)

   # plot nearest neighbour links
   b <- cells[m]
   arrows(cells$x, cells$y, b$x, b$y, angle=15, length=0.15, col="red")

   # find points which are the neighbour of their neighbour
   self <- (m[m] == seq(m))
   # plot them
   A <- cells[self]
   B <- cells[m[self]]
   plot(cells)
   segments(A$x, A$y, B$x, B$y)

   # nearest neighbours of each type
   head(nnwhich(ants, by=marks(ants)))
}
\author{Pavel Grabarnik
  \email{pavel.grabar@issp.serpukhov.su}
  and
  Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{spatial}
\keyword{math}

