\name{Extract.msr}
\alias{[.msr}
\title{Extract Subset of Signed or Vector Measure}
\description{
  Extract a subset of a signed measure or vector-valued measure.
}
\usage{
\method{[}{msr}(x, i, j, \dots)
}
\arguments{
  \item{x}{
    A signed or vector measure.
    An object of class \code{"msr"} (see \code{\link{msr}}).
  }
  \item{i}{
    Object defining the subregion or subset to be extracted.
    Either a spatial window (an object of class \code{"owin"}),
    or a pixel image with logical values,
    or any type of index that applies to a matrix.
  }
  \item{j}{
    Subset index selecting the vector coordinates to be extracted,
    if \code{x} is a vector-valued measure.
  }
  \item{\dots}{Ignored.}
} 
\value{
   An object of class \code{"msr"}.
}
\details{
  This operator extracts a subset of
  the data which determines the signed measure
  or vector-valued measure \code{x}. The result is another measure.
}
\seealso{
  \code{\link{msr}}
}
\examples{
   example(msr)
   rp[square(0.5)]
   rs[ , 2:3]
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au},
 Rolf Turner \email{r.turner@auckland.ac.nz}
 and Ege Rubak \email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{manip}
