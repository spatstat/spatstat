\name{pool.fv}
\alias{pool.fv}
\title{Pool Several Functions}
\description{
  Combine several summary functions into a single function.
}
\usage{
\method{pool}{fv}(..., weights=NULL, relabel=TRUE, variance=TRUE)
}
\arguments{
  \item{\dots}{
    Objects of class \code{"fv"}.
  }
  \item{weights}{
    Optional numeric vector of weights for the functions.
  }
  \item{relabel}{
    Logical value indicating whether the columns of the resulting function
    should be labelled to show that they were obtained by pooling.
  }
  \item{variance}{
    Logical value indicating whether to compute the sample variance
    and related terms.
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
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{htest}
\keyword{hplot}
\keyword{iteration}
