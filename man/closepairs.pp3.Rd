\name{closepairs.pp3}
\alias{closepairs.pp3}
\alias{crosspairs.pp3}
\title{
  Close Pairs of Points in 3 Dimensions
}
\description{
  Low-level functions to find all close pairs of points
  in three-dimensional point patterns.
}
\usage{
\method{closepairs}{pp3}(X, rmax, ordered=TRUE,
                         what=c("all", "indices"), \dots)

\method{crosspairs}{pp3}(X, Y, rmax, what=c("all", "indices"), \dots) 
}
\arguments{
  \item{X,Y}{
    Point patterns in three dimensions (objects of class \code{"pp3"}).
  }
  \item{rmax}{
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
    their \code{x,y,z} coordinates, and the distance between them.
    If \code{what="indices"} then only the indices \code{i,j} are returned.
  }
  \item{\dots}{Ignored.}
}
\details{
  These are the efficient low-level functions used by \pkg{spatstat}
  to find all close pairs of points in a three-dimensional point pattern
  or all close pairs between two point patterns in three dimensions.

  \code{closepairs(X,rmax)} identifies all pairs of neighbours 
  in the pattern \code{X} and returns them. The result is
  a list with the following components:
  \describe{
    \item{i}{Integer vector of indices of the first point in each pair.}
    \item{j}{Integer vector of indices of the second point in each pair.}
    \item{xi,yi,zi}{Coordinates of the first point in each pair.}
    \item{xj,yj,zj}{Coordinates of the second point in each pair.}
    \item{dx}{Equal to \code{xj-xi}}
    \item{dy}{Equal to \code{yj-yi}}
    \item{dz}{Equal to \code{zj-zi}}
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
  A list with components \code{i} and \code{j},
  and possibly other components as described under Details.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/},
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak.
}
\seealso{
  \code{\link{closepairs}}
}
\examples{
   X <- pp3(runif(10), runif(10), runif(10), box3(c(0,1)))
   Y <- pp3(runif(10), runif(10), runif(10), box3(c(0,1)))
   a <- closepairs(X, 0.1)
   b <- crosspairs(X, Y, 0.1)
}
\keyword{spatial}
\keyword{math}
