\name{FmultiInhom}
\alias{FmultiInhom}
\title{
  Inhomogeneous Marked F-Function
}
\description{
  For a marked point pattern, 
  estimate the inhomogeneous version of the multitype \eqn{F} function,
  effectively the cumulative distribution function of the distance from
  a fixed point to the nearest point in subset \eqn{J},
  adjusted for spatially varying intensity.
}
\usage{
  FmultiInhom(X, J,
              lambda = NULL, lambdaJ = NULL, lambdamin = NULL,
              \dots,
              r = NULL)
}
\arguments{
  \item{X}{
    A spatial point pattern (object of class \code{"ppp"}.
  }
  \item{J}{
    A subset index specifying the subset of points to which
    distances are measured. Any kind of subset index acceptable
    to \code{\link{[.ppp}}.
  }
  \item{lambda}{
    Intensity estimates for each point of \code{X}.
    A numeric vector of length equal to \code{npoints(X)}.
    Incompatible with \code{lambdaJ}.
  }
  \item{lambdaJ}{
    Intensity estimates for each point of \code{X[J]}.
    A numeric vector of length equal to \code{npoints(X[J])}.
    Incompatible with \code{lambda}.
  }
  \item{lambdamin}{
    A lower bound for the intensity,
    or at least a lower bound for the values in \code{lambdaJ}
    or \code{lambda[J]}.
  }
  \item{\dots}{
    Ignored. 
  }
  \item{r}{
    Vector of distance values at which the inhomogeneous \eqn{G}
    function should be estimated. There is a sensible default.
  }
}
\details{
  See Cronie and Van Lieshout (2015).
}
\value{
  Object of class \code{"fv"} containing the estimate of the
  inhomogeneous multitype \eqn{F} function.
}
\references{
  Cronie, O. and Van Lieshout, M.N.M. (2015)
  Summary statistics for inhomogeneous marked point processes.
  \emph{Annals of the Institute of Statistical Mathematics}
  DOI: 10.1007/s10463-015-0515-z
}
\author{
  Ottmar Cronie and Marie-Colette van Lieshout.
  Rewritten for \pkg{spatstat} by Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}.
}
\seealso{
  \code{\link{Finhom}}
}
\examples{
  X <- amacrine
  J <- (marks(X) == "off")
  mod <- ppm(X ~ marks * x)
  lam <- fitted(mod, dataonly=TRUE)
  lmin <- min(predict(mod)[["off"]]) * 0.9
  plot(FmultiInhom(X, J, lambda=lam, lambdamin=lmin))
}
\keyword{spatial}
\keyword{nonparametric}
