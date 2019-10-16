\name{Lcross}
\alias{Lcross}
\title{Multitype L-function (cross-type)}
\description{
  Calculates an estimate of the cross-type L-function
  for a multitype point pattern.
}
\usage{
  Lcross(X, i, j, ..., from, to, correction)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the cross-type \eqn{L} function
    \eqn{L_{ij}(r)}{Lij(r)} will be computed.
    It must be a multitype point pattern (a marked point pattern
    whose marks are a factor). See under Details.
  }
  \item{i}{The type (mark value)
    of the points in \code{X} from which distances are measured.
    A character string (or something that will be converted to a
    character string).
    Defaults to the first level of \code{marks(X)}.
  }
  \item{j}{The type (mark value)
    of the points in \code{X} to which distances are measured.
    A character string (or something that will be
    converted to a character string).
    Defaults to the second level of \code{marks(X)}.
  }
  \item{correction,\dots}{
    Arguments passed to \code{\link{Kcross}}.
  }
  \item{from,to}{
    An alternative way to specify \code{i} and \code{j} respectively.
  }
}
\details{
  The cross-type L-function is a transformation of the cross-type K-function,
  \deqn{L_{ij}(r) = \sqrt{\frac{K_{ij}(r)}{\pi}}}{Lij(r) = sqrt(Kij(r)/pi)}
  where \eqn{K_{ij}(r)}{Kij(r)} is the cross-type K-function
  from type \code{i} to type \code{j}.
  See \code{\link{Kcross}} for information
  about the cross-type K-function.

  The command \code{Lcross} first calls
  \code{\link{Kcross}} to compute the estimate of the cross-type K-function,
  and then applies the square root transformation.

  For a marked point pattern in which the points of type \code{i}
  are independent of the points of type \code{j},
  the theoretical value of the L-function is
  \eqn{L_{ij}(r) = r}{Lij(r) = r}.
  The square root also has the effect of stabilising
  the variance of the estimator, so that \eqn{L_{ij}}{Lij} is more appropriate
  for use in simulation envelopes and hypothesis tests.
}

\value{
  An object of class \code{"fv"}, see \code{\link{fv.object}},
  which can be plotted directly using \code{\link{plot.fv}}.

  Essentially a data frame containing columns
  \item{r}{the vector of values of the argument \eqn{r} 
    at which the function \eqn{L_{ij}}{Lij} has been  estimated
  }
  \item{theo}{the theoretical value \eqn{L_{ij}(r) = r}{Lij(r) = r}
    for a stationary Poisson process
  }
  together with columns named 
  \code{"border"}, \code{"bord.modif"},
  \code{"iso"} and/or \code{"trans"},
  according to the selected edge corrections. These columns contain
  estimates of the function \eqn{L_{ij}}{Lij} obtained by the edge corrections
  named.
}
\seealso{
  \code{\link{Kcross}},
  \code{\link{Ldot}},
  \code{\link{Lest}}
}
\examples{
 data(amacrine)
 L <- Lcross(amacrine, "off", "on")
 plot(L)
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{nonparametric}
