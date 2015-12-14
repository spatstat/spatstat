\name{K3est}
\Rdversion{1.1}
\alias{K3est}
\title{
  K-function of a Three-Dimensional Point Pattern
}
\description{
  Estimates the \eqn{K}-function from a three-dimensional point pattern.
}
\usage{
  K3est(X, \dots,
        rmax = NULL, nrval = 128,
        correction = c("translation", "isotropic"),
        ratio=FALSE)
}
\arguments{
  \item{X}{
    Three-dimensional point pattern (object of class \code{"pp3"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{rmax}{
    Optional. Maximum value of argument \eqn{r} for which
    \eqn{K_3(r)}{K3(r)} will be estimated. 
  }
  \item{nrval}{
    Optional. Number of values of \eqn{r} for which
    \eqn{K_3(r)}{K3(r)} will be estimated. A large value of \code{nrval}
    is required to avoid discretisation effects.
  }
  \item{correction}{
    Optional. Character vector specifying the edge correction(s)
    to be applied. See Details.
  }
  \item{ratio}{
    Logical. 
    If \code{TRUE}, the numerator and denominator of
    each edge-corrected estimate will also be saved,
    for use in analysing replicated point patterns.
  }
}
\details{
  For a stationary point process \eqn{\Phi}{Phi} in three-dimensional
  space, the three-dimensional \eqn{K} function
  is
  \deqn{
    K_3(r) = \frac 1 \lambda E(N(\Phi, x, r) \mid x \in \Phi)
  }{
    K3(r) = (1/lambda) E(N(Phi,x,r) | x in Phi)
  }
  where \eqn{\lambda}{lambda} is the intensity of the process
  (the expected number of points per unit volume) and
  \eqn{N(\Phi,x,r)}{N(Phi,x,r)} is the number of points of
  \eqn{\Phi}{Phi}, other than \eqn{x} itself, which fall within a
  distance \eqn{r} of \eqn{x}. This is the three-dimensional
  generalisation of Ripley's \eqn{K} function for two-dimensional
  point processes (Ripley, 1977).
  
  The three-dimensional point pattern \code{X} is assumed to be a
  partial realisation of a stationary point process \eqn{\Phi}{Phi}.
  The distance between each pair of distinct points is computed.
  The empirical cumulative distribution
  function of these values, with appropriate edge corrections, is
  renormalised to give the estimate of \eqn{K_3(r)}{K3(r)}.

  The available edge corrections are:
  \describe{
    \item{\code{"translation"}:}{
      the Ohser translation correction estimator
      (Ohser, 1983; Baddeley et al, 1993)
    }
    \item{\code{"isotropic"}:}{
      the three-dimensional counterpart of
      Ripley's isotropic edge correction (Ripley, 1977; Baddeley et al, 1993).
    }
  }
}
\value{
  A function value table (object of class \code{"fv"}) that can be
  plotted, printed or coerced to a data frame containing the function values.
}
\references{
  Baddeley, A.J, Moyeed, R.A., Howard, C.V. and Boyde, A. (1993)
  Analysis of a three-dimensional point pattern with replication.
  \emph{Applied Statistics} \bold{42}, 641--668.

  Ohser, J. (1983)
  On estimators for the reduced second moment measure of
  point processes. \emph{Mathematische Operationsforschung und
  Statistik, series Statistics}, \bold{14}, 63 -- 71.

  Ripley, B.D. (1977)
  Modelling spatial patterns (with discussion).
  \emph{Journal of the Royal Statistical Society, Series B},
  \bold{39}, 172 -- 212.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rana Moyeed.
}
\seealso{
  \code{\link{F3est}},
  \code{\link{G3est}},
  \code{\link{pcf3est}}
}
\examples{
  X <- rpoispp3(42)
  Z <- K3est(X)
  if(interactive()) plot(Z)
}
\keyword{spatial}
\keyword{nonparametric}
