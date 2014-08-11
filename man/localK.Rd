\name{localK}
\alias{localK}
\alias{localL}
\title{Neighbourhood density function}
\description{
  Computes the neighbourhood density function, a local version of
  the \eqn{K}-function or \eqn{L}-function,
  defined by Getis and Franklin (1987).
}
\usage{
  localK(X, ..., correction = "Ripley", verbose = TRUE, rvalue=NULL)
  localL(X, ..., correction = "Ripley", verbose = TRUE, rvalue=NULL)
}
\arguments{
  \item{X}{A point pattern (object of class \code{"ppp"}).}
  \item{\dots}{Ignored.}
  \item{correction}{String specifying the edge correction to be applied.
    Options are \code{"none"}, \code{"translate"}, \code{"translation"},
    \code{"Ripley"},
    \code{"isotropic"} or \code{"best"}.
    Only one correction may be specified.
  }
  \item{verbose}{Logical flag indicating whether to print progress
    reports during the calculation.
  }
  \item{rvalue}{Optional. A \emph{single} value of the distance argument
    \eqn{r} at which the function L or K should be computed.
  }
}
\details{
  The command \code{localL} computes the \emph{neighbourhood density function},
  a local version of the \eqn{L}-function (Besag's transformation of Ripley's
  \eqn{K}-function) that was proposed by Getis and Franklin (1987).
  The command \code{localK} computes the corresponding
  local analogue of the K-function.

  Given a spatial point pattern \code{X}, the neighbourhood density function
  \eqn{L_i(r)}{L[i](r)} associated with the \eqn{i}th point
  in \code{X} is computed by
  \deqn{
    L_i(r) = \sqrt{\frac a {(n-1) \pi} \sum_j e_{ij}}
  }{
    L[i](r) = sqrt( (a/((n-1)* pi)) * sum[j] e[i,j])
  }
  where the sum is over all points \eqn{j \neq i}{j != i} that lie
  within a distance \eqn{r} of the \eqn{i}th point, 
  \eqn{a} is the area of the observation window, \eqn{n} is the number
  of points in \code{X}, and \eqn{e_{ij}}{e[i,j]} is an edge correction
  term (as described in \code{\link{Kest}}).
  The value of \eqn{L_i(r)}{L[i](r)} can also be interpreted as one
  of the summands that contributes to the global estimate of the L
  function.

  By default, the function \eqn{L_i(r)}{L[i](r)} or
  \eqn{K_i(r)}{K[i](r)} is computed for a range of \eqn{r} values
  for each point \eqn{i}. The results are stored as a function value
  table (object of class \code{"fv"}) with a column of the table
  containing the function estimates for each point of the pattern
  \code{X}.

  Alternatively, if the argument \code{rvalue} is given, and it is a
  single number, then the function will only be computed for this value
  of \eqn{r}, and the results will be returned as a numeric vector,
  with one entry of the vector for each point of the pattern \code{X}.

  Inhomogeneous counterparts of \code{localK} and \code{localL}
  are computed by \code{localKinhom} and \code{localLinhom}.
}
\value{
  If \code{rvalue} is given, the result is a numeric vector
  of length equal to the number of points in the point pattern.

  If \code{rvalue} is absent, the result is 
  an object of class \code{"fv"}, see \code{\link{fv.object}},
  which can be plotted directly using \code{\link{plot.fv}}.
  Essentially a data frame containing columns
  \item{r}{the vector of values of the argument \eqn{r} 
    at which the function \eqn{K} has been  estimated
  }
  \item{theo}{the theoretical value \eqn{K(r) = \pi r^2}{K(r) = pi * r^2}
    or \eqn{L(r)=r} for a stationary Poisson process
  }
  together with columns containing the values of the
  neighbourhood density function for each point in the pattern.
  Column \code{i} corresponds to the \code{i}th point.
  The last two columns contain the \code{r} and \code{theo} values.
}
\references{
  Getis, A. and Franklin, J. (1987)
  Second-order neighbourhood analysis of mapped point patterns.
  \emph{Ecology} \bold{68}, 473--477.
}
\seealso{
  \code{\link{Kest}},
  \code{\link{Lest}},
  \code{\link{localKinhom}},
  \code{\link{localLinhom}}.
}
\examples{
  data(ponderosa)
  X <- ponderosa

  # compute all the local L functions
  L <- localL(X)

  # plot all the local L functions against r
  plot(L, main="local L functions for ponderosa", legend=FALSE)

  # plot only the local L function for point number 7
  plot(L, iso007 ~ r)
  
  # compute the values of L(r) for r = 12 metres
  L12 <- localL(X, rvalue=12)

  # Spatially interpolate the values of L12
  # Compare Figure 5(b) of Getis and Franklin (1987)
  X12 <- X \%mark\% L12
  Z <- Smooth(X12, sigma=5, dimyx=128)

  plot(Z, col=topo.colors(128), main="smoothed neighbourhood density")
  contour(Z, add=TRUE)
  points(X, pch=16, cex=0.5)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{nonparametric}

