\name{kaplan.meier}
\alias{kaplan.meier}
\title{Kaplan-Meier Estimator using Histogram Data}
\description{
  Compute the Kaplan-Meier estimator of a survival time distribution
  function, from histogram data
}
\usage{
  kaplan.meier(obs, nco, breaks, upperobs=0)
}
\arguments{
  \item{obs}{vector of \eqn{n} integers giving the histogram of
    all observations (censored or uncensored survival times)
  }
  \item{nco}{vector of \eqn{n} integers giving the histogram of
    uncensored observations (those survival times that are less than or
    equal to the censoring time)
  }
  \item{breaks}{Vector of \eqn{n+1} breakpoints which were used to form
    both histograms.
  }
  \item{upperobs}{
    Number of observations beyond the rightmost breakpoint, if any.
  }
}
\value{
  A list with two elements:
  \item{km}{Kaplan-Meier estimate of the survival time c.d.f. \eqn{F(t)}
  }
  \item{lambda}{corresponding Nelson-Aalen estimate of the
    hazard rate \eqn{\lambda(t)}{lambda(t)}
  }
  These are numeric vectors of length \eqn{n}.
}
\details{
  This function is needed mainly for internal use in \pkg{spatstat},
  but may be useful in other applications where you want to form the
  Kaplan-Meier estimator from a huge dataset.

  Suppose \eqn{T_i}{T[i]} are the survival times of individuals
  \eqn{i=1,\ldots,M} with unknown distribution function \eqn{F(t)}
  which we wish to estimate. Suppose these times are right-censored
  by random censoring times \eqn{C_i}{C[i]}.
  Thus the observations consist of right-censored survival times
  \eqn{\tilde T_i = \min(T_i,C_i)}{T*[i] = min(T[i],C[i])}
  and non-censoring indicators
  \eqn{D_i = 1\{T_i \le C_i\}}{D[i] = 1(T[i] <= C[i])}
  for each \eqn{i}.

  If the number of observations \eqn{M} is large, it is efficient to
  use histograms.
  Form the histogram \code{obs} of all observed times \eqn{\tilde T_i}{T*[i]}.
  That is, \code{obs[k]} counts the number of values 
  \eqn{\tilde T_i}{T*[i]} in the interval
  \code{(breaks[k],breaks[k+1]]} for \eqn{k > 1}
  and \code{[breaks[1],breaks[2]]} for \eqn{k = 1}.
  Also form the histogram \code{nco} of all uncensored times,
  i.e. those \eqn{\tilde T_i}{T*[i]} such that \eqn{D_i=1}{D[i]=1}.
  These two histograms are the arguments passed to \code{kaplan.meier}.
  
  The vectors \code{km} and \code{lambda} returned by \code{kaplan.meier}
  are (histogram approximations to) the Kaplan-Meier estimator
  of \eqn{F(t)} and its hazard rate \eqn{\lambda(t)}{lambda(t)}.
  Specifically, \code{km[k]} is an estimate of
  \code{F(breaks[k+1])}, and \code{lambda[k]} is an estimate of
  the average of \eqn{\lambda(t)}{lambda(t)} over the interval
  \code{(breaks[k],breaks[k+1])}.

  The histogram breaks must include \eqn{0}.
  If the histogram breaks do not span the range of the observations,
  it is important to count how many survival times
  \eqn{\tilde T_i}{T*[i]} exceed the rightmost breakpoint,
  and give this as the value \code{upperobs}.
}
\seealso{
  \code{\link{reduced.sample}},
  \code{\link{km.rs}}
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
  }
\keyword{spatial}
\keyword{nonparametric}

