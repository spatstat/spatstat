\name{as.data.frame.envelope}
\alias{as.data.frame.envelope}
\title{Coerce Envelope to Data Frame}
\description{
  Converts an envelope object to a data frame.
}
\usage{
\method{as.data.frame}{envelope}(x, \dots, simfuns=FALSE)
}
\arguments{
  \item{x}{Envelope object (class \code{"envelope"}).}
  \item{\dots}{Ignored.}
  \item{simfuns}{Logical value indicating whether the result should
    include the values of the simulated functions
    that were used to build the envelope. 
  }
}
\details{
  This is a method for the generic function \code{\link{as.data.frame}}
  for the class of envelopes (see \code{\link{envelope}}.

  The result is a data frame with columns
  containing the values of the function argument
  (usually named \code{r}), the function estimate for the original
  point pattern data (\code{obs}),
  the upper and lower envelope limits (\code{hi} and \code{lo}),
  and possibly additional columns.
  
  If \code{simfuns=TRUE}, the result also includes columns of values
  of the simulated functions that were used to compute the envelope.
  This is possible only when the envelope was computed with the
  argument \code{savefuns=TRUE} in the call to \code{\link{envelope}}.
}
\value{
  A data frame.
}
\examples{
  E <- envelope(cells, nsim=5, savefuns=TRUE)
  tail(as.data.frame(E))
  tail(as.data.frame(E, simfuns=TRUE))
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{manip}
