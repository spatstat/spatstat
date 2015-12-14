\name{pool.fv}
\alias{pool.fv}
\title{Pool Several Functions}
\description{
  Combine several summary functions into a single function.
}
\usage{
\method{pool}{fv}(..., weights=NULL)
}
\arguments{
  \item{\dots}{
    Objects of class \code{"fv"}.
  }
  \item{weights}{
    Optional numeric vector of weights for the functions.
  }
}
\details{
  The function \code{\link{pool}} is generic. This is the method for the
  class \code{"fv"} of summary functions. It is used to
  combine several estimates of the same function into a single function.

  Each of the arguments \code{\dots} must be an object of class
  \code{"fv"}. They must be compatible,
  in that they are estimates of the same function,
  and were computed using the same options.

  The sample mean and sample variance of the corresponding
  estimates will be computed.
}
\value{
  An object of class \code{"fv"}.
}
\seealso{
  \code{\link{pool}},
  \code{\link{pool.anylist}},
  \code{\link{pool.rat}}
}
\examples{
   K <- lapply(waterstriders, Kest, correction="iso")
   Kall <- pool(K[[1]], K[[2]], K[[3]])
   Kall <- pool(as.anylist(K))
   plot(Kall, cbind(pooliso, pooltheo) ~ r,
              shade=c("loiso", "hiiso"),
              main="Pooled K function of waterstriders")
}
\author{
Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}


Rolf Turner \email{r.turner@auckland.ac.nz}

and Ege Rubak \email{rubak@math.aau.dk}

}
\keyword{spatial}
\keyword{htest}
\keyword{hplot}
\keyword{iteration}
