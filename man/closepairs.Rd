\name{closepairs}
\alias{closepairs}
\alias{crosspairs}
\alias{closepairs.ppp}
\alias{crosspairs.ppp}
\alias{closepaircounts}
\alias{crosspaircounts}
\title{
  Close Pairs of Points
}
\description{
  Low-level functions to find all close pairs of points.
}
\usage{
closepaircounts(X, r)

crosspaircounts(X, Y, r)

closepairs(X, rmax, \dots)

\method{closepairs}{ppp}(X, rmax, ordered=TRUE,
                         what=c("all","indices"), \dots)

crosspairs(X, Y, rmax, \dots)

\method{crosspairs}{ppp}(X, Y, rmax, what=c("all", "indices"), \dots)
}
\arguments{
  \item{X,Y}{
    Point patterns (objects of class \code{"ppp"}).
  }
  \item{r,rmax}{
    Maximum distance between pairs of points to be counted as close pairs.
  }
  \item{ordered}{
    Logical value indicating whether all ordered pairs of close points
    should be returned. If \code{ordered=TRUE}, each pair will appear twice
    in the output, as \code{(i,j)} and again as \code{(j,i)}. If
    \code{ordered=FALSE}, then each pair will appear only once,
    as the pair \code{(i,j)} such that \code{i < j}.
  }
  \item{what}{
    String specifying the data to be returned for each close pair of points.
    If \code{what="all"} (the default) then the
    returned information includes the indices \code{i,j} of each pair,
    their \code{x,y} coordinates, and the distance between them.
    If \code{what="indices"} then only the indices \code{i,j} are returned.
  }
  \item{\dots}{Extra arguments, ignored by methods.}
}
\details{
  These are the efficient low-level functions used by \pkg{spatstat}
  to find all close pairs of points in a point pattern
  or all close pairs between two point patterns. 

  \code{closepaircounts(X,r)} counts the number of neighbours for
  each point in the pattern \code{X}. That is, for each point
  \code{X[i]}, it counts the number of other points \code{X[j]}
  with \code{j != i} such that \code{d(X[i],X[j]) <= r} where
  \code{d} denotes Euclidean distance. The result is an integer vector
  \code{v} such that \code{v[i]} is the number of neighbours of
  \code{X[i]}.

  \code{crosspaircounts(X,Y,r)} counts, for each point 
  in the pattern \code{X}, the number of neighbours in the pattern
  \code{Y}. That is, for each point
  \code{X[i]}, it counts the number of points \code{Y[j]}
  such that \code{d(X[i],X[j]) <= r}. The result is an integer vector
  \code{v} such that \code{v[i]} is the number of neighbours of
  \code{X[i]} in the pattern \code{Y}.

  \code{closepairs(X,rmax)} identifies all pairs of neighbours 
  in the pattern \code{X} and returns them. The result is
  a list with the following components:
  \describe{
    \item{i}{Integer vector of indices of the first point in each pair.}
    \item{j}{Integer vector of indices of the second point in each pair.}
    \item{xi,yi}{Coordinates of the first point in each pair.}
    \item{xj,yj}{Coordinates of the second point in each pair.}
    \item{dx}{Equal to \code{xj-xi}}
    \item{dy}{Equal to \code{yj-yi}}
    \item{d}{Euclidean distance between each pair of points.}
  }
  If \code{what="indices"} then only the components \code{i} and
  \code{j} are returned. This is slightly faster.

  \code{crosspairs(X,rmax)} identifies all pairs of neighbours
  \code{(X[i], Y[j])} between the patterns \code{X} and \code{Y},
  and returns them. The result is
  a list with the same format as for \code{closepairs}.
}
\section{Warning about accuracy}{
  The results of these functions may not agree exactly with
  the correct answer (as calculated by a human) and may not
  be consistent between different computers and different installations
  of \R. The discrepancies arise in marginal cases where the interpoint
  distance is equal to, or very close to, the threshold \code{rmax}.

  Floating-point numbers in a computer
  are not mathematical Real Numbers: they are approximations using
  finite-precision binary arithmetic.
  The approximation is accurate to a tolerance of about
  \code{.Machine$double.eps}.

  If the true interpoint distance \eqn{d} and the threshold \code{rmax}
  are equal, or if their difference is no more than \code{.Machine$double.eps},
  the result may be incorrect.
}
\value{
  For \code{closepaircounts} and \code{crosspaircounts}, an integer
  vector of length equal to the number of points in \code{X}.

  For \code{closepairs} and \code{crosspairs}, 
  a list with components \code{i} and \code{j},
  and possibly other components as described under Details.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{closepairs.pp3}} for the corresponding
  functions for 3D point patterns.
  
  \code{\link{Kest}}, \code{\link{Kcross}},
  \code{\link{nndist}}, \code{\link{nncross}},
  \code{\link{applynbd}}, \code{\link{markstat}}
  for functions which use these capabilities.
}
\examples{
   a <- closepaircounts(cells, 0.1)
   sum(a)

   Y <- split(amacrine)
   b <- crosspaircounts(Y$on, Y$off, 0.1)

   d <- closepairs(cells, 0.1)
   e <- crosspairs(Y$on, Y$off, 0.1)
}
\keyword{spatial}
\keyword{math}
