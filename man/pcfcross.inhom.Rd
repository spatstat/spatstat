\name{pcfcross.inhom}
\alias{pcfcross.inhom}
\title{
  Inhomogeneous Multitype Pair Correlation Function (Cross-Type)
}
\description{
  Estimates the inhomogeneous cross-type pair correlation function
  for a multitype point pattern.
}
\usage{
pcfcross.inhom(X, i, j, lambdaI = NULL, lambdaJ = NULL, ...,
               r = NULL, breaks = NULL,
               kernel="epanechnikov", bw=NULL, stoyan=0.15,
               correction = c("isotropic", "Ripley", "translate"),
               sigma = NULL, varcov = NULL)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the inhomogeneous
    cross-type pair correlation function
    \eqn{g_{ij}(r)}{g[i,j](r)}
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
  \item{j}{The type (mark value)
    of the points in \code{X} to which distances are measured.
    A character string (or something that will be
    converted to a character string).
    Defaults to the second level of \code{marks(X)}.
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
  \item{lambdaJ}{
    Optional.
    Values of the estimated intensity function of the points of type \code{j}.
    A numeric vector, pixel image or \code{function(x,y)}.
  }
  \item{r}{
    Vector of values for the argument \eqn{r} at which
    \eqn{g_{ij}(r)}{g[i,j](r)}
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
    \code{lambdaJ} is estimated by kernel smoothing.
  }
}
\details{
  The inhomogeneous cross-type pair correlation function
  \eqn{g_{ij}(r)}{g[i,j](r)}
  is a summary of the dependence between two types of points in a
  multitype spatial point process that does not have a uniform
  density of points.

  The best intuitive interpretation is the following: the probability
  \eqn{p(r)} of finding two points, of types \eqn{i} and \eqn{j}
  respectively, at locations \eqn{x} and \eqn{y}
  separated by a distance \eqn{r} is equal to
  \deqn{
    p(r) = \lambda_i(x) lambda_j(y) g(r) \,{\rm d}x \, {\rm d}y
  }{
    p(r) = lambda[i](x) * lambda[j](y) * g(r) dx dy
  }
  where \eqn{\lambda_i}{lambda[i]} is the intensity function
  of the process of points of type \eqn{i}.
  For a multitype Poisson point process,
  this probability is 
  \eqn{p(r) = \lambda_i(x) \lambda_j(y)}{p(r) = lambda[i](x) * lambda[j](y)}
  so   \eqn{g_{ij}(r) = 1}{g[i,j](r) = 1}.

  The command \code{pcfcross.inhom} estimates the inhomogeneous
  pair correlation using a modified version of
  the algorithm in \code{\link{pcf.ppp}}.

  If the arguments \code{lambdaI} and \code{lambdaJ} are missing or
  null, they are estimated from \code{X} by kernel smoothing using a
  leave-one-out estimator. 
}
\value{
  A function value table (object of class \code{"fv"}).
  Essentially a data frame containing the variables
  \item{r}{
    the vector of values of the argument \eqn{r} 
    at which the inhomogeneous cross-type pair correlation function
    \eqn{g_{ij}(r)}{g[i,j](r)} has been  estimated
  }
  \item{theo}{vector of values equal to 1,
    the theoretical value of \eqn{g_{ij}(r)}{g[i,j](r)}
    for the Poisson process
  }
  \item{trans}{vector of values of \eqn{g_{ij}(r)}{g[i,j](r)}
    estimated by translation correction
  }
  \item{iso}{vector of values of \eqn{g_{ij}(r)}{g[i,j](r)}
    estimated by Ripley isotropic correction
  }
  as required.
}
\seealso{
  \code{\link{pcf.ppp}}, 
  \code{\link{pcfinhom}}, 
  \code{\link{pcfcross}},
  \code{\link{pcfdot.inhom}}
}
\examples{
  data(amacrine)
  plot(pcfcross.inhom(amacrine, "on", "off", stoyan=0.1),
       legendpos="bottom")
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{nonparametric}
