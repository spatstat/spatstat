\name{envelopeArray}
\alias{envelopeArray}
\title{
  Array of Simulation Envelopes of Summary Function
}
\description{
  Compute an array of simulation envelopes using
  a summary function that returns an array of curves.
}
\usage{
envelopeArray(X, fun, \dots, dataname = NULL, verb = FALSE, reuse = TRUE)
}
\arguments{
  \item{X}{
    Object containing point pattern data.
    A point pattern (object of class
    \code{"ppp"}, \code{"lpp"}, \code{"pp3"} or \code{"ppx"})
    or a fitted point process model
    (object of class \code{"ppm"}, \code{"kppm"} or \code{"lppm"}).
  }
  \item{fun}{
    Function that computes the desired summary statistic
    for a point pattern. The result of \code{fun} should be a
    function array (object of class \code{"fasp"}).
  }
  \item{\dots}{
  Arguments passed to \code{\link{envelope}} to control the simulations,
  or passed to \code{fun} when evaluating the function.
  }
  \item{dataname}{
    Optional character string name for the data.
  }
  \item{verb}{
    Logical value indicating whether to print progress reports.
  }
  \item{reuse}{
    Logical value indicating whether the envelopes in each panel
    should be based on the same set of simulated patterns
    (\code{reuse=TRUE}, the default)
    or on different, independent sets of simulated
    patterns (\code{reuse=FALSE}).
  }
}
\details{
  This command is the counterpart of \code{\link{envelope}}
  when the function \code{fun} that is evaluated on each simulated point pattern
  will return an object of class \code{"fasp"} representing an array of
  summary functions.

  Simulated point patterns are generated according to the
  rules described for \code{\link{envelope}}. In brief, 
  if \code{X} is a point pattern, the algorithm generates
  simulated point patterns of the same kind, according to complete
  spatial randomness. If \code{X} is a fitted model, the algorithm
  generates simulated point patterns according to this model.

  For each simulated point pattern \code{Y}, the function \code{fun}
  is invoked. The result \code{Z <- fun(Y, ...)} should be an object of
  class \code{"fasp"} representing an array of summary functions.
  The dimensions of the array \code{Z} should be the same
  for each simulated pattern \code{Y}.

  This algorithm finds the simulation envelope of the summary functions
  in each cell of the array. 
}
\value{
   An object of class \code{"fasp"} representing
   an array of envelopes.
}
\author{
   \spatstatAuthors.
}
\seealso{
  \code{\link{envelope}}, \code{\link{alltypes}}.
}
\examples{
  A <- envelopeArray(finpines, markcrosscorr, nsim=9)
  plot(A)
}
\keyword{spatial}
\keyword{nonparametric}
\keyword{iteration}
