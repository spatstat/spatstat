\name{pcfdot.inhom}
\alias{pcfdot.inhom}
\title{
  Inhomogeneous Multitype Pair Correlation Function (Type-i-To-Any-Type)
}
\description{
  Estimates the inhomogeneous multitype pair correlation function
  (from type \eqn{i} to any type)
  for a multitype point pattern.
}
\usage{
pcfdot.inhom(X, i, lambdaI = NULL, lambdadot = NULL, ...,
               r = NULL, breaks = NULL,
               kernel="epanechnikov", bw=NULL, stoyan=0.15,
               correction = c("isotropic", "Ripley", "translate"),
               sigma = NULL, varcov = NULL)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the inhomogeneous
    multitype pair correlation function
    \eqn{g_{i\bullet}(r)}{g[i.](r)}
    will be computed.
    It must be a multitype point pattern (a marked point pattern
    whose marks are a factor). 
  }
  \item{i}{The type (mark value)
    of the points in \code{X} from which distances are measured.
    A character string (or something that will be converted to a
    character string).
    Defaults to the first level of \code{marks(X)}.
  }
  \item{lambdaI}{
    Optional.
    Values of the estimated intensity function of the points of type \code{i}.
    Either a vector giving the intensity values
    at the points of type \code{i},
    a pixel image (object of class \code{"im"}) giving the
    intensity values at all locations, or a \code{function(x,y)} which
    can be evaluated to give the intensity value at any location.
  }
  \item{lambdadot}{
    Optional.
    Values of the estimated intensity function of the point pattern \code{X}.
    A numeric vector, pixel image or \code{function(x,y)}.
  }
  \item{r}{
    Vector of values for the argument \eqn{r} at which
    \eqn{g_{i\bullet}(r)}{g[i.](r)}
    should be evaluated. There is a sensible default.
  }
  \item{breaks}{
	This argument is for internal use only.
  }
  \item{kernel}{
    Choice of smoothing kernel, passed to \code{\link{density.default}}.
  }
  \item{bw}{
    Bandwidth for smoothing kernel, passed to \code{\link{density.default}}.
  }
  \item{\dots}{
    Other arguments passed to the kernel density estimation 
    function \code{\link{density.default}}.
  }
  \item{stoyan}{
    Bandwidth coefficient; see Details.
  }
  \item{correction}{
    Choice of edge correction.
  }
  \item{sigma,varcov}{
    Optional arguments passed to  \code{\link{density.ppp}}
    to control the smoothing bandwidth, when \code{lambdaI} or
    \code{lambdadot} is estimated by kernel smoothing.
  }
}
\details{
  The inhomogeneous multitype (type \eqn{i} to any type)
  pair correlation function
  \eqn{g_{i\bullet}(r)}{g[i.](r)}
  is a summary of the dependence between different types of points in a
  multitype spatial point process that does not have a uniform
  density of points.

  The best intuitive interpretation is the following: the probability
  \eqn{p(r)} of finding a point of type \eqn{i} at location \eqn{x}
  and another point of any type at location \eqn{y},
  where \eqn{x} and \eqn{y} are separated by a distance \eqn{r},
  is equal to
  \deqn{
    p(r) = \lambda_i(x) lambda(y) g(r) \,{\rm d}x \, {\rm d}y
  }{
    p(r) = lambda[i](x) * lambda(y) * g(r) dx dy
  }
  where \eqn{\lambda_i}{lambda[i]} is the intensity function
  of the process of points of type \eqn{i}, and 
  where \eqn{\lambda}{lambda} is the intensity function of the points
  of all types.
  For a multitype Poisson point process, this probability is 
  \eqn{p(r) = \lambda_i(x) \lambda(y)}{p(r) = lambda[i](x) * lambda(y)}
  so   \eqn{g_{i\bullet}(r) = 1}{g[i.](r) = 1}.

  The command \code{pcfdot.inhom} estimates the inhomogeneous
  multitype pair correlation using a modified version of
  the algorithm in \code{\link{pcf.ppp}}.

  If the arguments \code{lambdaI} and \code{lambdadot} are missing or
  null, they are estimated from \code{X} by kernel smoothing using a
  leave-one-out estimator. 
}
\value{
  A function value table (object of class \code{"fv"}).
  Essentially a data frame containing the variables
  \item{r}{
    the vector of values of the argument \eqn{r} 
    at which the inhomogeneous multitype pair correlation function
    \eqn{g_{i\bullet}(r)}{g[i.](r)} has been  estimated
  }
  \item{theo}{vector of values equal to 1,
    the theoretical value of \eqn{g_{i\bullet}(r)}{g[i.](r)}
    for the Poisson process
  }
  \item{trans}{vector of values of \eqn{g_{i\bullet}(r)}{g[i.](r)}
    estimated by translation correction
  }
  \item{iso}{vector of values of \eqn{g_{i\bullet}(r)}{g[i.](r)}
    estimated by Ripley isotropic correction
  }
  as required.
}
\seealso{
  \code{\link{pcf.ppp}}, 
  \code{\link{pcfinhom}}, 
  \code{\link{pcfdot}},
  \code{\link{pcfcross.inhom}}
}
\examples{
  data(amacrine)
  plot(pcfdot.inhom(amacrine, "on", stoyan=0.1), legendpos="bottom")
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{nonparametric}
