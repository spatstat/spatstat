\name{km.rs}
\alias{km.rs}
\title{Kaplan-Meier and Reduced Sample Estimator using Histograms}
\description{
  Compute the Kaplan-Meier and Reduced Sample estimators of a
  survival time distribution function, using histogram techniques
}
\usage{
  km.rs(o, cc, d, breaks)
}
\arguments{
  \item{o}{vector of observed survival times
  }
  \item{cc}{vector of censoring times
  }
  \item{d}{vector of non-censoring indicators
  }
  \item{breaks}{Vector of breakpoints to be used to form histograms.
  }
}
\value{
  A list with five elements
  \item{rs}{Reduced-sample estimate of the survival time c.d.f. \eqn{F(t)}
  }
  \item{km}{Kaplan-Meier estimate of the survival time c.d.f. \eqn{F(t)}
  }
  \item{hazard}{corresponding Nelson-Aalen estimate of the
    hazard rate \eqn{\lambda(t)}{lambda(t)}
  }
  \item{r}{values of \eqn{t} for which \eqn{F(t)} is estimated
  }
  \item{breaks}{the breakpoints vector
  }
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

  The arguments to this function are 
  vectors \code{o}, \code{cc}, \code{d}
  of observed values of \eqn{\tilde T_i}{T*[i]}, \eqn{C_i}{C[i]}
  and \eqn{D_i}{D[i]} respectively.
  The function computes histograms and forms the reduced-sample
  and Kaplan-Meier estimates of  \eqn{F(t)} by
  invoking the functions \code{\link{kaplan.meier}}
  and \code{\link{reduced.sample}}.
  This is efficient if the lengths of \code{o}, \code{cc}, \code{d}
  (i.e. the number of observations) is large.

  The vectors \code{km} and \code{hazard} returned by \code{kaplan.meier}
  are (histogram approximations to) the Kaplan-Meier estimator
  of \eqn{F(t)} and its hazard rate \eqn{\lambda(t)}{lambda(t)}.
  Specifically, \code{km[k]} is an estimate of
  \code{F(breaks[k+1])}, and \code{lambda[k]} is an estimate of
  the average of \eqn{\lambda(t)}{lambda(t)} over the interval
  \code{(breaks[k],breaks[k+1])}. This approximation is exact only if the
  survival times are discrete and the 
  histogram breaks are fine enough to ensure that each interval
  \code{(breaks[k],breaks[k+1])} contains only one possible value of
  the survival time. 

  The vector \code{rs} is the reduced-sample estimator,
  \code{rs[k]} being the reduced sample estimate of \code{F(breaks[k+1])}.
  This value is exact, i.e. the use of histograms does not introduce any
  approximation error in the reduced-sample estimator.
}
\seealso{
  \code{\link{reduced.sample}},
  \code{\link{kaplan.meier}}
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
  }
\keyword{spatial}
\keyword{nonparametric}

