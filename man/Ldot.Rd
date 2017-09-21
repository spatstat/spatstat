\name{Ldot}
\alias{Ldot}
\title{Multitype L-function (i-to-any)}
\description{
  Calculates an estimate of the multitype L-function
  (from type \code{i} to any type)
  for a multitype point pattern.
}
\usage{
  Ldot(X, i, ..., from)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the dot-type \eqn{L} function
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
  \item{\dots}{
    Arguments passed to \code{\link{Kdot}}.
  }
  \item{from}{An alternative way to specify \code{i}.}
}
\details{
  This command computes 
  \deqn{L_{i\bullet}(r) = \sqrt{\frac{K_{i\bullet}(r)}{\pi}}}{Li.(r) = sqrt(Ki.(r)/pi)}
  where \eqn{K_{i\bullet}(r)}{Ki.(r)} is the multitype \eqn{K}-function
  from points of type \code{i} to points of any type.
  See \code{\link{Kdot}} for information
  about \eqn{K_{i\bullet}(r)}{Ki.(r)}.

  The command \code{Ldot} first calls
  \code{\link{Kdot}} to compute the estimate of the \code{i}-to-any
  \eqn{K}-function, and then applies the square root transformation.

  For a marked Poisson point process,
  the theoretical value of the L-function is
  \eqn{L_{i\bullet}(r) = r}{Li.(r) = r}.
  The square root also has the effect of stabilising
  the variance of the estimator, so that \eqn{L_{i\bullet}}{Li.}
  is more appropriate
  for use in simulation envelopes and hypothesis tests.
}

\value{
  An object of class \code{"fv"}, see \code{\link{fv.object}},
  which can be plotted directly using \code{\link{plot.fv}}.

  Essentially a data frame containing columns
  \item{r}{the vector of values of the argument \eqn{r} 
    at which the function \eqn{L_{i\bullet}}{Li.} has been  estimated
  }
  \item{theo}{the theoretical value \eqn{L_{i\bullet}(r) = r}{Li.(r) = r}
    for a stationary Poisson process
  }
  together with columns named 
  \code{"border"}, \code{"bord.modif"},
  \code{"iso"} and/or \code{"trans"},
  according to the selected edge corrections. These columns contain
  estimates of the function \eqn{L_{i\bullet}}{Li.}
  obtained by the edge corrections named.
}
\seealso{
  \code{\link{Kdot}},
  \code{\link{Lcross}},
  \code{\link{Lest}}
}
\examples{
 data(amacrine)
 L <- Ldot(amacrine, "off")
 plot(L)
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{nonparametric}
