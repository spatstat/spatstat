\name{localKcross.inhom}
\alias{localKcross.inhom}
\alias{localLcross.inhom}
\title{Inhomogeneous Multitype K Function}
\description{
  Computes spatially-weighted versions of the
  the local multitype \eqn{K}-function or \eqn{L}-function.
}
\usage{
  localKcross.inhom(X, from, to,
              lambdaFrom=NULL, lambdaTo=NULL,
              \dots, rmax = NULL,
              correction = "Ripley", sigma=NULL, varcov=NULL,
              lambdaX=NULL, update=TRUE, leaveoneout=TRUE)
  localLcross.inhom(X, from, to,
              lambdaFrom=NULL, lambdaTo=NULL, \dots, rmax = NULL)
}
\arguments{
  \item{X}{A point pattern (object of class \code{"ppp"}).}
  \item{from}{
    Type of points from which distances should be measured.
    A single value;
    one of the possible levels of \code{marks(X)},
    or an integer indicating which level.
  }
  \item{to}{
    Type of points to which distances should be measured.
    A single value; 
    one of the possible levels of \code{marks(X)},
    or an integer indicating which level.
  }
  \item{lambdaFrom,lambdaTo}{
    Optional.
    Values of the estimated intensity function
    for the points of type \code{from} and \code{to}, respectively.
    Each argument should be either a vector giving the intensity values
    at the required points,
    a pixel image (object of class \code{"im"}) giving the
    intensity values at all locations, a fitted point process model
    (object of class \code{"ppm"}) or a \code{function(x,y)} which
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
  \item{sigma, varcov}{
    Optional arguments passed to \code{\link{density.ppp}} to control
    the kernel smoothing procedure for estimating \code{lambdaFrom}
    and \code{lambdaTo}, if they are missing.
  }
  \item{lambdaX}{
    Optional. 
    Values of the estimated intensity function
    for all points of \code{X}.
    Either a vector giving the intensity values
    at each point of \code{X},
    a pixel image (object of class \code{"im"}) giving the
    intensity values at all locations, a fitted point process model
    (object of class \code{"ppm"}) or a \code{function(x,y)} which
    can be evaluated to give the intensity value at any location.
  }
  \item{update}{
    Logical value indicating what to do when \code{lambdaFrom},
    \code{lambdaTo} or \code{lambdaX} is a fitted model
    (class \code{"ppm"}, \code{"kppm"} or \code{"dppm"}).
    If \code{update=TRUE} (the default),
    the model will first be refitted to the data \code{X}
    (using \code{\link{update.ppm}} or \code{\link{update.kppm}})
    before the fitted intensity is computed.
    If \code{update=FALSE}, the fitted intensity of the
    model will be computed without re-fitting it to \code{X}.
  }
  \item{leaveoneout}{
    Logical value (passed to \code{\link{density.ppp}} or
    \code{\link{fitted.ppm}}) specifying whether to use a
    leave-one-out rule when calculating the intensity.
  }
}
\details{
  The functions \code{localKcross.inhom} and \code{localLcross.inhom}
  are inhomogeneous or weighted versions of the
  local multitype \eqn{K} and \eqn{L} functions implemented in
  \code{\link{localKcross}} and \code{\link{localLcross}}.

  Given a multitype spatial point pattern \code{X},
  and two designated types \code{from} and \code{to},
  the local multitype \eqn{K} function is
  defined for each point \code{X[i]} that belongs to type \code{from},
  and is computed by
  \deqn{
    K_i(r) = \sqrt{\frac 1 \pi \sum_j \frac{e_{ij}}{\lambda_j}}
  }{
    K[i](r) = sqrt( (1/pi) * sum[j] e[i,j]/lambda[j])
  }
  where the sum is over all points \eqn{j \neq i}{j != i}
  of type \code{to} that lie
  within a distance \eqn{r} of the \eqn{i}th point, 
  \eqn{\lambda_j}{\lambda[j]} is the estimated intensity of the
  point pattern at the point \eqn{j},
  and \eqn{e_{ij}}{e[i,j]} is an edge correction
  term (as described in \code{\link{Kest}}).

  The function 
  \eqn{K_i(r)}{K[i](r)} is computed for a range of \eqn{r} values
  for each point \eqn{i}. The results are stored as a function value
  table (object of class \code{"fv"}) with a column of the table
  containing the function estimates for each point of the pattern
  \code{X} of type \code{from}.

  The corresponding \eqn{L} function 
  \eqn{L_i(r)}{L[i](r)} is computed by applying the
  transformation 
  \eqn{L(r) = \sqrt{K(r)/(2\pi)}}{L(r) = sqrt(K(r)/(2*pi))}.
}
\value{
  An object of class \code{"fv"}, see \code{\link{fv.object}},
  which can be plotted directly using \code{\link{plot.fv}}.
  Essentially a data frame containing columns
  \item{r}{the vector of values of the argument \eqn{r} 
    at which the function \eqn{K} has been  estimated
  }
  \item{theo}{the theoretical value \eqn{K(r) = \pi r^2}{K(r) = pi * r^2}
    or \eqn{L(r)=r} for a stationary Poisson process
  }
  together with columns containing the values of the
  neighbourhood density function for each point in the pattern
  of type \code{from}.
  The last two columns contain the \code{r} and \code{theo} values.
}
\seealso{
  \code{\link{Kinhom}},
  \code{\link{Linhom}},
  \code{\link{localK}},
  \code{\link{localL}}.
}
\examples{
  X <- amacrine

  # compute all the local L functions
  L <- localLcross.inhom(X)

  # plot all the local L functions against r
  plot(L, main="local L functions for ponderosa", legend=FALSE)

  # plot only the local L function for point number 7
  plot(L, iso007 ~ r)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{nonparametric}

