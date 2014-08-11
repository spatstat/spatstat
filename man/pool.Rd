\name{pool}
\alias{pool}
\title{
  Pool Data 
}
\description{
  Pool the data from several objects of the same class.
}
\usage{
pool(...)
}
\arguments{
  \item{\dots}{
    Objects of the same type.
  }
}
\details{
  The function \code{pool} is generic. There are methods for several
  classes, listed below.
  
  \code{pool} is used to combine the data from several objects of the same type,
  and to compute statistics based on the combined dataset.
  It may be used to pool the estimates obtained from replicated datasets.
  It may also be used in high-performance computing applications,
  when the objects \code{\dots} have been computed on different processors
  or in different batch runs, and we wish to combine them.
}
\value{
  An object of the same class as the arguments \code{\dots}.
}
\seealso{
  \code{\link{pool.envelope}},
  \code{\link{pool.fasp}},
  \code{\link{pool.rat}}
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
