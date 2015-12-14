\name{reduced.sample}
\alias{reduced.sample}
\title{Reduced Sample Estimator using Histogram Data}
\description{
  Compute the Reduced Sample estimator of a survival time distribution
  function, from histogram data
}
\usage{
  reduced.sample(nco, cen, ncc, show=FALSE, uppercen=0)
}
\arguments{
  \item{nco}{vector of counts giving the histogram of
    uncensored observations (those survival times that are less than or
    equal to the censoring time)
  }
  \item{cen}{vector of counts giving the histogram of
    censoring times
  }
  \item{ncc}{vector of counts giving the histogram of
    censoring times for the uncensored observations only 
  }
  \item{uppercen}{
    number of censoring times greater than the rightmost
    histogram breakpoint (if there are any)
  }
  \item{show}{Logical value controlling the amount of detail
    returned by the function value (see below)
  }
}
\value{
  If \code{show = FALSE}, a numeric vector giving the values of
  the reduced sample estimator.
  If \code{show=TRUE}, a list with three components which are
  vectors of equal length,
  \item{rs}{Reduced sample estimate of the survival time c.d.f. \eqn{F(t)}
  }
  \item{numerator}{numerator of the reduced sample estimator
  }
  \item{denominator}{denominator of the reduced sample estimator
  }
}
\details{
  This function is needed mainly for internal use in \pkg{spatstat},
  but may be useful in other applications where you want to form the
  reduced sample estimator from a huge dataset.

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
  Form the histogram \code{cen} of all censoring times \eqn{C_i}{C[i]}.
  That is, \code{obs[k]} counts the number of values 
  \eqn{C_i}{C[i]} in the interval
  \code{(breaks[k],breaks[k+1]]} for \eqn{k > 1}
  and \code{[breaks[1],breaks[2]]} for \eqn{k = 1}.
  Also form the histogram \code{nco} of all uncensored times,
  i.e. those \eqn{\tilde T_i}{T*[i]} such that \eqn{D_i=1}{D[i]=1},
  and the histogram of all censoring times for which the survival time
  is uncensored,
  i.e. those \eqn{C_i}{C[i]} such that \eqn{D_i=1}{D[i]=1}.
  These three histograms are the arguments passed to \code{kaplan.meier}.

  The return value \code{rs} is the reduced-sample estimator
  of the distribution function \eqn{F(t)}. Specifically,
  \code{rs[k]} is the reduced sample estimate of \code{F(breaks[k+1])}.
  The value is exact, i.e. the use of histograms does not introduce any
  approximation error.

  Note that, for the results to be valid, either the histogram breaks
  must span the censoring times, or the number of censoring times
  that do not fall in a histogram cell must have been counted in
  \code{uppercen}.
}
\seealso{
  \code{\link{kaplan.meier}},
  \code{\link{km.rs}}
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
  }
\keyword{spatial}
\keyword{nonparametric}

