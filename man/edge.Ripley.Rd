\name{edge.Ripley}
\alias{edge.Ripley}
\alias{rmax.Ripley}
\title{
  Ripley's Isotropic Edge Correction 
}
\description{
  Computes Ripley's isotropic edge correction weights
  for a point pattern.
}
\usage{
edge.Ripley(X, r, W = Window(X), method = c("C", "interpreted"),
            maxweight = 100, internal=list())

rmax.Ripley(W)
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"}).
  }
  \item{W}{
    Window for which the edge correction is required.
  }
  \item{r}{
    Vector or matrix of interpoint distances for which the edge correction
    should be computed.
  }
  \item{method}{
    Choice of algorithm. Either \code{"interpreted"} or \code{"C"}.
    This is needed only for debugging purposes.
  }
  \item{maxweight}{
    Maximum permitted value of the edge correction weight.
  }
  \item{internal}{For developer use only.}
}
\details{
  The function \code{edge.Ripley}
  computes Ripley's (1977) isotropic edge correction
  weight, which is used in estimating the \eqn{K} function and in many
  other contexts.

  The function \code{rmax.Ripley} computes the maximum value of
  distance \eqn{r} for which the isotropic edge correction
  estimate of \eqn{K(r)} is valid.
    
  For a single point \eqn{x} in a window \eqn{W},
  and a distance \eqn{r > 0}, the isotropic edge correction weight
  is
  \deqn{
    e(u, r) = \frac{2\pi r}{\mbox{length}(c(u,r) \cap W)}
  }{
    e(u, r) = 2 * \pi * r/length(intersection(c(u,r), W))
  }
  where \eqn{c(u,r)} is the circle of radius \eqn{r} centred at the
  point \eqn{u}. The denominator is the length of the overlap between
  this circle and the window \eqn{W}.

  The function \code{edge.Ripley} computes this edge correction weight
  for each point in the point pattern \code{X} and for each
  corresponding distance value in the vector or matrix \code{r}.
  
  If \code{r} is a vector, with one entry for each point in
  \code{X}, then the result is a vector containing the
  edge correction weights \code{e(X[i], r[i])} for each \code{i}.

  If \code{r} is a matrix, with one row for each point in \code{X},
  then the result is a matrix whose \code{i,j} entry gives the
  edge correction weight \code{e(X[i], r[i,j])}.
  For example \code{edge.Ripley(X, pairdist(X))} computes all the
  edge corrections required for the \eqn{K}-function.

  If any value of the edge correction weight exceeds \code{maxwt},
  it is set to \code{maxwt}.

  The function \code{rmax.Ripley} computes the smallest distance \eqn{r}
  such that it is possible to draw a circle of radius \eqn{r}, centred
  at a point of \code{W}, such that the circle does not intersect the
  interior of \code{W}. 
}
\value{
  A numeric vector or matrix.
}
\references{
  Ripley, B.D. (1977)
  Modelling spatial patterns (with discussion).
  \emph{Journal of the Royal Statistical Society, Series B},
  \bold{39}, 172 -- 212.
}
\seealso{
  \code{\link{edge.Trans}},
  \code{\link{rmax.Trans}},
  \code{\link{Kest}}
}
\examples{
  v <- edge.Ripley(cells, pairdist(cells))

  rmax.Ripley(Window(cells))
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{nonparametric}

