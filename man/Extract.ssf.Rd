\name{[.ssf}
\alias{[.ssf}
\title{
  Subset of spatially sampled function
}
\description{
  Extract a subset of the data for a spatially sampled function.
}
\usage{
 \method{[}{ssf}(x, i, j, ..., drop)
}
\arguments{
  \item{x}{
    Object of class \code{"ssf"}.
  }  
  \item{i}{
    Subset index applying to the locations where the function is sampled.
  }
  \item{j}{
   Subset index applying to the columns (variables) measured at each location.
  }
  \item{\dots, drop}{
   Ignored.
 }
}
\details{
   This is the subset operator for the class \code{"ssf"}.
}
\value{
  Another object of class \code{"ssf"}.
}
\author{
\adrian.
}
\seealso{
\code{\link{ssf}}, 
\code{\link{with.ssf}}
}
\examples{
  f <- ssf(cells, data.frame(d=nndist(cells), i=1:42))
  f
  f[1:10,]
  f[ ,1]
}
\keyword{spatial}
\keyword{manip}
