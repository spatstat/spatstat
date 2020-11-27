\name{closepaircounts}
\alias{closepaircounts}
\alias{crosspaircounts}
\title{
  Count Close Pairs of Points
}
\description{
  Low-level functions to count the number of close pairs of points.
}
\usage{
closepaircounts(X, r)

crosspaircounts(X, Y, r)
}
\arguments{
  \item{X,Y}{
    Point patterns (objects of class \code{"ppp"}).
  }
  \item{r}{
    Maximum distance between pairs of points to be counted as close pairs.
  }
}
\details{
  These are the efficient low-level functions used by \pkg{spatstat}
  to count close pairs of points in a point pattern
  or between two point patterns. 

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
  An integer vector of length equal to the number of points in \code{X}.
}
\author{
  \adrian
  and \rolf
}
\seealso{
  \code{\link{closepairs}} to identify all close pairs of points.
}
\examples{
   a <- closepaircounts(cells, 0.1)
   sum(a)
   Y <- split(amacrine)
   b <- crosspaircounts(Y$on, Y$off, 0.1)
}
\keyword{spatial}
\keyword{math}
