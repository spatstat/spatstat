\name{GmultiInhom}
\alias{GmultiInhom}
\title{
  Inhomogeneous Marked G-Function
}
\description{
  For a marked point pattern, 
  estimate the inhomogeneous version of the multitype \eqn{G} function,
  effectively the cumulative distribution function of the distance from
  a point in subset \eqn{I} to the nearest point in subset \eqn{J},
  adjusted for spatially varying intensity.
}
\usage{
  GmultiInhom(X, I, J,
             lambda = NULL, lambdaI = NULL, lambdaJ = NULL,
             lambdamin = NULL, \dots,
             r = NULL,
             ReferenceMeasureMarkSetI = NULL,
             ratio = FALSE)
}
\arguments{
  \item{X}{
    A spatial point pattern (object of class \code{"ppp"}.
  }
  \item{I}{
    A subset index specifying the subset of points \emph{from} which
    distances are measured. Any kind of subset index acceptable
    to \code{\link{[.ppp}}.
  }
  \item{J}{
    A subset index specifying the subset of points \emph{to} which
    distances are measured. Any kind of subset index acceptable
    to \code{\link{[.ppp}}.
  }
  \item{lambda}{
    Intensity estimates for each point of \code{X}.
    A numeric vector of length equal to \code{npoints(X)}.
    Incompatible with \code{lambdaI,lambdaJ}.
  }
  \item{lambdaI}{
    Intensity estimates for each point of \code{X[I]}.
    A numeric vector of length equal to \code{npoints(X[I])}.
    Incompatible with \code{lambda}.
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
  \item{ReferenceMeasureMarkSetI}{
    Optional. The total measure of the mark set. A positive number.
  }
  \item{ratio}{
    Logical value indicating whether to save ratio information.
  }
}
\details{
  See Cronie and Van Lieshout (2015).
}
\value{
  Object of class \code{"fv"} containing the estimate of the
  inhomogeneous multitype \eqn{G} function.
}
\references{
  Cronie, O. and Van Lieshout, M.N.M. (2015)
  Summary statistics for inhomogeneous marked point processes.
  \emph{Annals of the Institute of Statistical Mathematics}
  DOI: 10.1007/s10463-015-0515-z
}
\author{
  Ottmar Cronie and Marie-Colette van Lieshout.
  Rewritten for \pkg{spatstat} by \adrian.
}
\seealso{
  \code{\link{Ginhom}},
  \code{\link{Gmulti}}
}
\examples{
  X <- amacrine
  I <- (marks(X) == "on")
  J <- (marks(X) == "off")
  mod <- ppm(X ~ marks * x)
  lam <- fitted(mod, dataonly=TRUE)
  lmin <- min(predict(mod)[["off"]]) * 0.9
  plot(GmultiInhom(X, I, J, lambda=lam, lambdamin=lmin))
}
\keyword{spatial}
\keyword{nonparametric}
