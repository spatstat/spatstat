\name{Ginhom}
\alias{Ginhom}
\title{
  Inhomogeneous Nearest Neighbour Function
}
\description{
  Estimates the inhomogeneous nearest neighbour function \eqn{G} of
  a non-stationary point pattern.
}
\usage{
  Ginhom(X, lambda = NULL, lmin = NULL, ...,
        sigma = NULL, varcov = NULL,
        r = NULL, breaks = NULL, ratio = FALSE)
}
\arguments{
  \item{X}{
    The observed data point pattern,
    from which an estimate of the inhomogeneous \eqn{G} function
    will be computed.
    An object of class \code{"ppp"}
    or in a format recognised by \code{\link{as.ppp}()}
  }
  \item{lambda}{
    Optional.
    Values of the estimated intensity function.
    Either a vector giving the intensity values
    at the points of the pattern \code{X},
    a pixel image (object of class \code{"im"}) giving the
    intensity values at all locations, a fitted point process model
    (object of class \code{"ppm"}) or a \code{function(x,y)} which
    can be evaluated to give the intensity value at any location.
  }
  \item{lmin}{
    Optional. The minimum possible value of the intensity
    over the spatial domain. A positive numerical value.
  }
  \item{sigma,varcov}{
    Optional arguments passed to  \code{\link{density.ppp}}
    to control the smoothing bandwidth, when \code{lambda} is
    estimated by kernel smoothing.
  }
  \item{\dots}{
    Extra arguments passed to \code{\link{as.mask}} to control
    the pixel resolution, or passed to \code{\link{density.ppp}}
    to control the smoothing bandwidth.
}
  \item{r}{
    vector of values for the argument \eqn{r} at which
    the inhomogeneous \eqn{K} function
    should be evaluated. Not normally given by the user;
    there is a sensible default.
  }
  \item{breaks}{
    This argument is for internal use only.
  }
  \item{ratio}{
    Logical. 
    If \code{TRUE}, the numerator and denominator of
    the estimate will also be saved,
    for use in analysing replicated point patterns.
  }
}
\details{
  This command computes estimates of the 
  inhomogeneous \eqn{G}-function (van Lieshout, 2010)
  of a point pattern. It is the counterpart, for inhomogeneous
  spatial point patterns, of the nearest-neighbour distance
  distribution function \eqn{G} 
  for homogeneous point patterns computed by \code{\link{Gest}}.

  The argument \code{X} should be a point pattern
  (object of class \code{"ppp"}).

  The inhomogeneous \eqn{G} function is computed
  using the border correction, equation (7) in Van Lieshout (2010).
  
  The argument \code{lambda} should supply the
  (estimated) values of the intensity function \eqn{\lambda}{lambda}
  of the point process. It may be either
  \describe{
    \item{a numeric vector}{
      containing the values
      of the intensity function at the points of the pattern \code{X}.
    }
    \item{a pixel image}{
      (object of class \code{"im"})
      assumed to contain the values of the intensity function
      at all locations in the window. 
    }
    \item{a fitted point process model}{
      (object of class \code{"ppm"})
      whose fitted \emph{trend} can be used as the fitted intensity.
    }
    \item{a function}{
      which can be evaluated to give values of the intensity at
      any locations.
    }
    \item{omitted:}{
      if \code{lambda} is omitted, then it will be estimated using
      a `leave-one-out' kernel smoother.
    }
  }
  If \code{lambda} is a numeric vector, then its length should
  be equal to the number of points in the pattern \code{X}.
  The value \code{lambda[i]} is assumed to be the 
  the (estimated) value of the intensity
  \eqn{\lambda(x_i)}{lambda(x[i])} for
  the point \eqn{x_i}{x[i]} of the pattern \eqn{X}.
  Each value must be a positive number; \code{NA}'s are not allowed.

  If \code{lambda} is a pixel image, the domain of the image should
  cover the entire window of the point pattern. If it does not (which
  may occur near the boundary because of discretisation error),
  then the missing pixel values 
  will be obtained by applying a Gaussian blur to \code{lambda} using
  \code{\link{blur}}, then looking up the values of this blurred image
  for the missing locations. 
  (A warning will be issued in this case.)

  If \code{lambda} is a function, then it will be evaluated in the
  form \code{lambda(x,y)} where \code{x} and \code{y} are vectors
  of coordinates of the points of \code{X}. It should return a numeric
  vector with length equal to the number of points in \code{X}.

  If \code{lambda} is omitted, then it will be estimated using
  a `leave-one-out' kernel smoother,
  as described in Baddeley, \ifelse{latex}{\out{M\o ller}}{Moller}
  and Waagepetersen (2000).  The estimate \code{lambda[i]} for the
  point \code{X[i]} is computed by removing \code{X[i]} from the
  point pattern, applying kernel smoothing to the remaining points using
  \code{\link{density.ppp}}, and evaluating the smoothed intensity
  at the point \code{X[i]}. The smoothing kernel bandwidth is controlled
  by the arguments \code{sigma} and \code{varcov}, which are passed to
  \code{\link{density.ppp}} along with any extra arguments.
}
\value{
  An object of class \code{"fv"}, see \code{\link{fv.object}},
  which can be plotted directly using \code{\link{plot.fv}}.
}
\references{
  Baddeley, A.,
  \ifelse{latex}{\out{M\o ller}}{Moller}, J. and Waagepetersen, R. (2000)
  Non- and semiparametric estimation of interaction in
  inhomogeneous point patterns.
  \emph{Statistica Neerlandica} \bold{54}, 329--350.

  Van Lieshout, M.N.M. and Baddeley, A.J. (1996)
  A nonparametric measure of spatial interaction in point patterns.
  \emph{Statistica Neerlandica} \bold{50}, 344--361.

  Van Lieshout, M.N.M. (2010)
  A J-function for inhomogeneous point processes.
  \emph{Statistica Neerlandica} \bold{65}, 183--201.
}
\seealso{
  \code{\link{Finhom}},
  \code{\link{Jinhom}},
  \code{\link{Gest}}
}
\examples{
  \dontrun{
    plot(Ginhom(swedishpines, sigma=bw.diggle, adjust=2))
  }
  plot(Ginhom(swedishpines, sigma=10))
}
\author{
  Original code by Marie-Colette van Lieshout.
  C implementation and R adaptation by Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{spatial}
\keyword{nonparametric}
