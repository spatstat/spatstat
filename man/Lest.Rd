\name{Lest}
\alias{Lest}
\title{L-function}
\description{
  Calculates an estimate of the \eqn{L}-function (Besag's
  transformation of Ripley's \eqn{K}-function)
  for a spatial point pattern.
}
\usage{
  Lest(X, ...)
}
\arguments{
  \item{X}{
    The observed point pattern, 
    from which an estimate of \eqn{L(r)} will be computed.
    An object of class \code{"ppp"}, or data
    in any format acceptable to \code{\link{as.ppp}()}.
  }
  \item{\dots}{
    Other arguments passed to \code{\link{Kest}}
    to control the estimation procedure.
  }
}
\details{
  This command computes an estimate of the \eqn{L}-function
  for the spatial point pattern \code{X}.
  The \eqn{L}-function is a transformation of Ripley's \eqn{K}-function,
  \deqn{L(r) = \sqrt{\frac{K(r)}{\pi}}}{L(r) = sqrt(K(r)/pi)}
  where \eqn{K(r)} is the \eqn{K}-function.

  See \code{\link{Kest}} for information
  about Ripley's \eqn{K}-function. The transformation to \eqn{L} was
  proposed by Besag (1977).

  The command \code{Lest} first calls
  \code{\link{Kest}} to compute the estimate of the \eqn{K}-function,
  and then applies the square root transformation.

  For a completely random (uniform Poisson) point pattern,
  the theoretical value of the \eqn{L}-function is \eqn{L(r) = r}.
  The square root also has the effect of stabilising
  the variance of the estimator, so that \eqn{L(r)} is more appropriate
  for use in simulation envelopes and hypothesis tests.

  See \code{\link{Kest}} for the list of arguments.
}
\section{Variance approximations}{
  If the argument \code{var.approx=TRUE} is given, the return value
  includes columns \code{rip} and \code{ls} containing approximations
  to the variance of \eqn{\hat L(r)}{Lest(r)} under CSR.
  These are obtained by the delta method from the variance
  approximations described in \code{\link{Kest}}.
}

\value{
  An object of class \code{"fv"}, see \code{\link{fv.object}},
  which can be plotted directly using \code{\link{plot.fv}}.

  Essentially a data frame containing columns
  \item{r}{the vector of values of the argument \eqn{r} 
    at which the function \eqn{L} has been  estimated
  }
  \item{theo}{the theoretical value \eqn{L(r) = r}
    for a stationary Poisson process
  }
  together with columns named 
  \code{"border"}, \code{"bord.modif"},
  \code{"iso"} and/or \code{"trans"},
  according to the selected edge corrections. These columns contain
  estimates of the function \eqn{L(r)} obtained by the edge corrections
  named.
}
\references{
  Besag, J. (1977) 
  Discussion of Dr Ripley's paper.
  \emph{Journal of the Royal Statistical Society, Series B},
  \bold{39}, 193--195.
}
\seealso{
  \code{\link{Kest}},
  \code{\link{pcf}}
}
\examples{
 data(cells)
 L <- Lest(cells)
 plot(L, main="L function for cells")
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{nonparametric}
