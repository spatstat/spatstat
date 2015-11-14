\name{pool.fasp}
\alias{pool.fasp}
\title{
  Pool Data from Several Function Arrays
}
\description{
  Pool the simulation data from several function arrays
  (objects of class \code{"fasp"})
  and compute a new function array.
}
\usage{
\method{pool}{fasp}(...)
}
\arguments{
  \item{\dots}{
    Objects of class \code{"fasp"}.
  }
}
\details{
  The function \code{\link{pool}} is generic. This is the method for the
  class \code{"fasp"} of function arrays. It is used to
  combine the simulation data from several arrays of simulation envelopes
  and to compute a new array of envelopes based on the combined data.

  Each of the arguments \code{\dots} must be a function array
  (object of class \code{"fasp"}) containing simulation envelopes.
  This is typically created by running the command
  \code{\link{alltypes}} with the arguments
  \code{envelope=TRUE} and \code{savefuns=TRUE}.
  This ensures that each object is an array of simulation envelopes,
  and that each envelope contains the simulated data
  (summary function values) that were used to construct the envelope.

  The simulated data are extracted from each object and combined.
  A new array of envelopes is computed from the combined set of simulations.

  Warnings or errors will be issued if the objects \code{\dots}
  appear to be incompatible. However, the code is not smart enough to
  decide whether it is sensible to pool the data.
}
\value{
  An object of class \code{"fasp"}.
}
\seealso{
  \code{\link{fasp}},
  \code{\link{alltypes}},
  \code{\link{pool.envelope}},
  \code{\link{pool}}
}
\examples{
   data(amacrine)
   A1 <- alltypes(amacrine,"K",nsim=9,envelope=TRUE,savefuns=TRUE)
   A2 <- alltypes(amacrine,"K",nsim=10,envelope=TRUE,savefuns=TRUE)
   pool(A1, A2)
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{htest}
\keyword{hplot}
\keyword{iteration}
