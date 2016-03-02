\name{edge.Trans}
\alias{edge.Trans}
\alias{rmax.Trans}
\title{
  Translation Edge Correction
}
\description{
  Computes Ohser and Stoyan's translation edge correction weights
  for a point pattern.
}
\usage{
edge.Trans(X, Y = X, W = Window(X),
      exact = FALSE, paired = FALSE,
      ..., 
      trim = spatstat.options("maxedgewt"),
      dx=NULL, dy=NULL,
      give.rmax=FALSE, gW=NULL)

rmax.Trans(W, g=setcov(W))
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
  \item{\dots}{Ignored.}
  \item{trim}{
    Maximum permitted value of the edge correction weight.
  }
  \item{dx,dy}{
    Alternative data giving the \eqn{x} and \eqn{y} coordinates
    of the vector differences between the points.
    Incompatible with \code{X} and \code{Y}. See Details.
  }
  \item{give.rmax}{
    Logical. If \code{TRUE}, also compute the value of
    \code{rmax.Trans(W)} and return it as an attribute
    of the result.
  }
  \item{g, gW}{
    Optional. Set covariance of \code{W}, if it has already been
    computed. Not required if \code{W} is a rectangle.
  }
}
\details{
  The function \code{edge.Trans}
  computes Ohser and Stoyan's translation edge correction
  weight, which is used in estimating the \eqn{K} function and in many
  other contexts.

  The function \code{rmax.Trans} computes the maximum value of
  distance \eqn{r} for which the translation edge correction
  estimate of \eqn{K(r)} is valid.
  
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

  The arguments \code{dx} and \code{dy} can be provided as
  an alternative to \code{X} and \code{Y}.
  If \code{paired=TRUE} then \code{dx,dy} should be vectors of equal length
  such that the vector difference of the \eqn{i}th pair is
  \code{c(dx[i], dy[i])}. If \code{paired=FALSE} then 
  \code{dx,dy} should be matrices of the same dimensions,
  such that the vector difference between \code{X[i]} and \code{Y[j]} is
  \code{c(dx[i,j], dy[i,j])}. The argument \code{W} is needed.

  The value of \code{rmax.Trans} is the shortest distance from the
  origin \eqn{(0,0)} to the boundary of the support of
  the set covariance function of \code{W}. It is computed by pixel
  approximation using \code{\link{setcov}}, unless \code{W} is a
  rectangle, when \code{rmax.Trans(W)} is the length of the
  shortest side of the rectangle.
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
  \code{\link{rmax.Trans}},
  \code{\link{edge.Ripley}},
  \code{\link{setcov}},
  \code{\link{Kest}}
}
\examples{
  v <- edge.Trans(cells)
  rmax.Trans(Window(cells))
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  and Rolf Turner \email{r.turner@auckland.ac.nz}.
}
\keyword{spatial}
\keyword{nonparametric}

