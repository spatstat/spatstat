\name{edge.Trans}
\alias{edge.Trans}
\title{
  Translation Edge Correction
}
\description{
  Computes Ohser and Stoyan's translation edge correction weights
  for a point pattern.
}
\usage{
edge.Trans(X, Y = X, W = X$window,
      exact = FALSE, paired = FALSE,
      trim = spatstat.options("maxedgewt"))
}
\arguments{
  \item{X,Y}{
    Point patterns (objects of class \code{"ppp"}).
  }
  \item{W}{
    Window for which the edge correction is required.
  }
  \item{exact}{
    Logical. If \code{TRUE}, a slow algorithm will be used
    to compute the exact value. If \code{FALSE}, a fast algorithm
    will be used to compute the approximate value.
  }
  \item{paired}{
    Logical value indicating whether \code{X} and \code{Y}
    are paired. If \code{TRUE}, compute
    the edge correction for corresponding points
    \code{X[i], Y[i]} for all \code{i}.
    If \code{FALSE}, compute the edge correction for
    each possible pair of points \code{X[i], Y[j]}
    for all \code{i} and \code{j}.
  }
  \item{trim}{
    Maximum permitted value of the edge correction weight.
  }
}
\details{
  This function computes Ohser and Stoyan's translation edge correction
  weight, which is used in estimating the \eqn{K} function and in many
  other contexts.

  For a pair of points \eqn{x} and \eqn{y} in a window \eqn{W},
  the translation edge correction weight
  is
  \deqn{
    e(u, r) = \frac{\mbox{area}(W)}{\mbox{area}(W \cap (W + y - x))}
  }{
    e(u, r) = area(W) / area(intersect(W, W + y - x))
  }
  where \eqn{W + y - x} is the result of shifting the window \eqn{W}
  by the vector \eqn{y - x}. The denominator is the area of the overlap between
  this shifted window and the original window.

  The function \code{edge.Trans} computes this edge correction weight.
  If \code{paired=TRUE}, then \code{X} and \code{Y} should contain the
  same number of points. The result is a vector containing the
  edge correction weights \code{e(X[i], Y[i])} for each \code{i}.

  If \code{paired=FALSE}, 
  then the result is a matrix whose \code{i,j} entry gives the
  edge correction weight \code{e(X[i], Y[j])}.

  Computation is exact if the window is a rectangle.
  Otherwise,
  \itemize{
    \item if \code{exact=TRUE}, the edge
    correction weights are computed exactly using 
    \code{\link{overlap.owin}}, which can be quite slow.
    \item if \code{exact=FALSE} (the default),
    the weights are computed rapidly by evaluating the
    set covariance function \code{\link{setcov}}
    using the Fast Fourier Transform.
  }
  If any value of the edge correction weight exceeds \code{trim},
  it is set to \code{trim}. 
}
\value{
  Numeric vector or matrix.
}
\references{
  Ohser, J. (1983)
  On estimators for the reduced second moment measure of
  point processes. \emph{Mathematische Operationsforschung und
  Statistik, series Statistics}, \bold{14}, 63 -- 71.
}
\seealso{
  \code{\link{edge.Ripley}},
  \code{\link{setcov}},
  \code{\link{Kest}}
}
\examples{
  v <- edge.Trans(cells)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{nonparametric}

