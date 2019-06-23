\name{localKinhom}
\alias{localKinhom}
\alias{localLinhom}
\title{Inhomogeneous Neighbourhood Density Function}
\description{
  Computes spatially-weighted versions of the
  the local \eqn{K}-function or \eqn{L}-function.
}
\usage{
  localKinhom(X, lambda, ..., rmax = NULL,
              correction = "Ripley", verbose = TRUE, rvalue=NULL,
              sigma = NULL, varcov = NULL, update=TRUE, leaveoneout=TRUE)
  localLinhom(X, lambda, ..., rmax = NULL, 
              correction = "Ripley", verbose = TRUE, rvalue=NULL,
              sigma = NULL, varcov = NULL, update=TRUE, leaveoneout=TRUE)
}
\arguments{
  \item{X}{A point pattern (object of class \code{"ppp"}).}
  \item{lambda}{
    Optional.
    Values of the estimated intensity function.
    Either a vector giving the intensity values
    at the points of the pattern \code{X},
    a pixel image (object of class \code{"im"}) giving the
    intensity values at all locations, a fitted point process model
    (object of class \code{"ppm"} or \code{"kppm"} or \code{"dppm"})
    or a \code{function(x,y)} which
    can be evaluated to give the intensity value at any location.
  }
  \item{\dots}{
    Extra arguments. Ignored if \code{lambda} is present.
    Passed to \code{\link{density.ppp}} if \code{lambda} is omitted.
  }
  \item{rmax}{
    Optional. Maximum desired value of the argument \eqn{r}.
  }
  \item{correction}{
    String specifying the edge correction to be applied.
    Options are \code{"none"}, \code{"translate"}, \code{"Ripley"},
    \code{"translation"}, \code{"isotropic"} or \code{"best"}.
    Only one correction may be specified.
  }
  \item{verbose}{Logical flag indicating whether to print progress
    reports during the calculation.
  }
  \item{rvalue}{Optional. A \emph{single} value of the distance argument
    \eqn{r} at which the function L or K should be computed.
  }
  \item{sigma, varcov}{
    Optional arguments passed to \code{\link{density.ppp}} to control
    the kernel smoothing procedure for estimating \code{lambda},
    if \code{lambda} is missing.
  }
  \item{leaveoneout}{
    Logical value (passed to \code{\link{density.ppp}} or
    \code{\link{fitted.ppm}}) specifying whether to use a
    leave-one-out rule when calculating the intensity.
  }
  \item{update}{
    Logical value indicating what to do when \code{lambda} is a fitted model
    (class \code{"ppm"}, \code{"kppm"} or \code{"dppm"}).
    If \code{update=TRUE} (the default),
    the model will first be refitted to the data \code{X}
    (using \code{\link{update.ppm}} or \code{\link{update.kppm}})
    before the fitted intensity is computed.
    If \code{update=FALSE}, the fitted intensity of the
    model will be computed without re-fitting it to \code{X}.
  }
}
\details{
  The functions \code{localKinhom} and \code{localLinhom}
  are inhomogeneous or weighted versions of the
  neighbourhood density function implemented in
  \code{\link{localK}} and \code{\link{localL}}.

  Given a spatial point pattern \code{X}, the
  inhomogeneous neighbourhood density function
  \eqn{L_i(r)}{L[i](r)} associated with the \eqn{i}th point
  in \code{X} is computed by
  \deqn{
    L_i(r) = \sqrt{\frac 1 \pi \sum_j \frac{e_{ij}}{\lambda_j}}
  }{
    L[i](r) = sqrt( (1/pi) * sum[j] e[i,j]/lambda[j])
  }
  where the sum is over all points \eqn{j \neq i}{j != i} that lie
  within a distance \eqn{r} of the \eqn{i}th point, 
  \eqn{\lambda_j}{\lambda[j]} is the estimated intensity of the
  point pattern at the point \eqn{j},
  and \eqn{e_{ij}}{e[i,j]} is an edge correction
  term (as described in \code{\link{Kest}}).
  The value of \eqn{L_i(r)}{L[i](r)} can also be interpreted as one
  of the summands that contributes to the global estimate of the
  inhomogeneous L function (see \code{\link{Linhom}}).

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
\seealso{
  \code{\link{Kinhom}},
  \code{\link{Linhom}},
  \code{\link{localK}},
  \code{\link{localL}}.
}
\examples{
  data(ponderosa)
  X <- ponderosa

  # compute all the local L functions
  L <- localLinhom(X)

  # plot all the local L functions against r
  plot(L, main="local L functions for ponderosa", legend=FALSE)

  # plot only the local L function for point number 7
  plot(L, iso007 ~ r)
  
  # compute the values of L(r) for r = 12 metres
  L12 <- localL(X, rvalue=12)
}
\author{
  Mike Kuhn,
  \adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{nonparametric}

