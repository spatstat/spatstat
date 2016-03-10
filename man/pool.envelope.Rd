\name{pool.envelope}
\alias{pool.envelope}
\title{
  Pool Data from Several Envelopes
}
\description{
  Pool the simulation data from several simulation envelopes
  (objects of class \code{"envelope"})
  and compute a new envelope.
}
\usage{
\method{pool}{envelope}(..., savefuns=FALSE, savepatterns=FALSE)
}
\arguments{
  \item{\dots}{
    Objects of class \code{"envelope"}.
  }
  \item{savefuns}{
    Logical flag indicating whether to save all the simulated
    function values.
  }
  \item{savepatterns}{
    Logical flag indicating whether to save all the simulated
    point patterns.
  }
}
\details{
  The function \code{\link{pool}} is generic. This is the method for the
  class \code{"envelope"} of simulation envelopes. It is used to
  combine the simulation data from several simulation envelopes
  and to compute an envelope based on the combined data.

  Each of the arguments \code{\dots} must be an object of class
  \code{"envelope"}. These envelopes must be compatible,
  in that they are envelopes for the same function,
  and were computed using the same options.

  \itemize{
    \item In normal use, each envelope object will have been 
    created by running the command \code{\link{envelope}}
    with the argument \code{savefuns=TRUE}.
    This ensures that each object contains the simulated data
    (summary function values for the simulated point patterns)
    that were used to construct the envelope.

    The simulated data are extracted from each object and combined.
    A new envelope is computed from the combined set of simulations.
    \item
    Alternatively, if each envelope object was created
    by running \code{\link{envelope}} with \code{VARIANCE=TRUE},
    then the saved functions are not required.
    
    The sample means and sample variances from each envelope
    will be pooled. A new envelope is computed from the pooled
    mean and variance.
  }

  Warnings or errors will be issued if the envelope objects \code{\dots}
  appear to be incompatible. Apart from these basic checks,
  the code is not smart enough to decide whether it is sensible
  to pool the data.

  To modify the envelope parameters or the type of envelope that is
  computed, first pool the envelope data using \code{pool.envelope},
  then use \code{\link{envelope.envelope}} to modify the envelope
  parameters.
}
\value{
  An object of class \code{"envelope"}.
}
\seealso{
  \code{\link{envelope}},
  \code{\link{envelope.envelope}},
  \code{\link{pool}},
  \code{\link{pool.fasp}}
}
\examples{
   E1 <- envelope(cells, Kest, nsim=10, savefuns=TRUE)
   E2 <- envelope(cells, Kest, nsim=20, savefuns=TRUE)
   pool(E1, E2)

   V1 <- envelope(E1, VARIANCE=TRUE)
   V2 <- envelope(E2, VARIANCE=TRUE)
   pool(V1, V2)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{htest}
\keyword{hplot}
\keyword{iteration}
