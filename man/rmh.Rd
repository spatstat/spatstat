\name{rmh}
\alias{rmh}
\title{Simulate point patterns using the Metropolis-Hastings algorithm.}
\description{
  Generic function for running the Metropolis-Hastings algorithm
  to produce simulated realisations of a point process model.
}

\usage{rmh(model, \dots)}

\arguments{
  \item{model}{The point process model to be simulated.
  }
  \item{\dots}{Further arguments controlling the simulation.
  }
}

\details{
  The Metropolis-Hastings algorithm can be used to
  generate simulated realisations from a wide range of
  spatial point processes. For caveats, see below.
  
  The function \code{rmh} is generic; it has methods
  \code{\link{rmh.ppm}} (for objects of class \code{"ppm"})
  and  \code{\link{rmh.default}} (the default).
  The actual implementation of the Metropolis-Hastings algorithm is
  contained in \code{\link{rmh.default}}.
  For details of its use, see 
  \code{\link{rmh.ppm}} or \code{\link{rmh.default}}.

  [If the model is a Poisson process, then Metropolis-Hastings
  is not used; the Poisson model is generated directly
  using \code{\link{rpoispp}} or \code{\link{rmpoispp}}.]

  In brief, the Metropolis-Hastings algorithm is a Markov Chain,
  whose states are spatial point patterns, and whose limiting
  distribution is the desired point process. After
  running the algorithm for a very large number of iterations,
  we may regard the state of the algorithm as a realisation
  from the desired point process.

  However, there are difficulties in deciding whether the
  algorithm has run for ``long enough''. The convergence of the
  algorithm may indeed be extremely slow. No guarantees of
  convergence are given!

  While it is fashionable to decry the Metropolis-Hastings algorithm
  for its poor convergence and other properties, it has the advantage
  of being easy to implement for a wide range of models.
}

\section{Warning}{
As of version 1.22-1 of \code{spatstat} a subtle change was
made to \code{rmh.default()}.  We had noticed that the results
produced were sometimes not ``scalable'' in that two models,
differing in effect only by the units in which distances are
measured and starting from the same seed, gave different results.
This was traced to an idiosyncracy of floating point arithmetic.
The code of \code{rmh.default()} has been changed so that the
results produced by \code{rmh} are now scalable.  The downside of
this is that code which users previously ran may now give results
which are different from what they formerly were.

In order to recover former behaviour (so that previous results
can be reproduced) set \code{spatstat.options(scalable=FALSE)}.
See the last example in the help for \code{\link{rmh.default}}.
}

\value{
  A point pattern, in the form of an object of class \code{"ppp"}.
  See \code{\link{rmh.default}} for details.
}

\seealso{
  \code{\link{rmh.default}}
}

\examples{
    # See examples in rmh.default and rmh.ppm
}

\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{datagen}

