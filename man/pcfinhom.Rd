\name{pcfinhom}
\alias{pcfinhom}
\title{
  Inhomogeneous Pair Correlation Function
}
\description{
  Estimates the inhomogeneous pair correlation function of
  a point pattern using kernel methods.
}
\usage{
pcfinhom(X, lambda = NULL, ..., r = NULL,
         kernel = "epanechnikov", bw = NULL, stoyan = 0.15,
         correction = c("translate", "Ripley"),
         divisor = c("r", "d"),
         renormalise = TRUE, normpower=1,
         reciplambda = NULL,
         sigma = NULL, varcov = NULL)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
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
  \item{r}{
    Vector of values for the argument \eqn{r} at which \eqn{g(r)} 
    should be evaluated. There is a sensible default.
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
    Bandwidth coefficient; see \code{\link{pcf.ppp}}.
  }
  \item{correction}{
    Choice of edge correction.
  }
  \item{divisor}{
    Choice of divisor in the estimation formula:
    either \code{"r"} (the default) or \code{"d"}.
    See \code{\link{pcf.ppp}}.
  }
  \item{renormalise}{
    Logical. Whether to renormalise the estimate. See Details.
  }
  \item{normpower}{
    Integer (usually either 1 or 2).
    Normalisation power. See Details.
  }
  \item{reciplambda}{
    Alternative to \code{lambda}.
    Values of the estimated \emph{reciprocal} \eqn{1/\lambda}{1/lambda}
    of the intensity function.
    Either a vector giving the reciprocal intensity values
    at the points of the pattern \code{X},
    a pixel image (object of class \code{"im"}) giving the
    reciprocal intensity values at all locations,
    or a \code{function(x,y)} which can be evaluated to give the
    reciprocal intensity value at any location.
  }
  \item{sigma,varcov}{
    Optional arguments passed to  \code{\link{density.ppp}}
    to control the smoothing bandwidth, when \code{lambda} is
    estimated by kernel smoothing.
  }
}
\details{
  The inhomogeneous pair correlation function \eqn{g_{\rm inhom}(r)}{ginhom(r)}
  is a summary of the dependence between points in a spatial point
  process that does not have a uniform density of points.

  The best intuitive interpretation is the following: the probability
  \eqn{p(r)} of finding two points at locations \eqn{x} and \eqn{y}
  separated by a distance \eqn{r} is equal to
  \deqn{
    p(r) = \lambda(x) lambda(y) g(r) \,{\rm d}x \, {\rm d}y
  }{
    p(r) = lambda(x) * lambda(y) * g(r) dx dy
  }
  where \eqn{\lambda}{lambda} is the intensity function
  of the point process.
  For a Poisson point process with intensity function
  \eqn{\lambda}{lambda}, this probability is 
  \eqn{p(r) = \lambda(x) \lambda(y)}{p(r) = lambda(x) * lambda(y)}
  so \eqn{g_{\rm inhom}(r) = 1}{ginhom(r) = 1}.

  The inhomogeneous pair correlation function 
  is related to the inhomogeneous \eqn{K} function through
  \deqn{
    g_{\rm inhom}(r) = \frac{K'_{\rm inhom}(r)}{2\pi r}
  }{
    ginhom(r) = Kinhom'(r)/ ( 2 * pi * r) 
  }
  where \eqn{K'_{\rm inhom}(r)}{Kinhom'(r)}
  is the derivative of \eqn{K_{\rm inhom}(r)}{Kinhom(r)}, the
  inhomogeneous \eqn{K} function. See \code{\link{Kinhom}} for information
  about \eqn{K_{\rm inhom}(r)}{Kinhom(r)}.

  The command \code{pcfinhom} estimates the inhomogeneous
  pair correlation using a modified version of
  the algorithm in \code{\link{pcf.ppp}}.   
  
  If \code{renormalise=TRUE} (the default), then the estimates 
  are multiplied by \eqn{c^{\mbox{normpower}}}{c^normpower} where 
  \eqn{
    c = \mbox{area}(W)/\sum (1/\lambda(x_i)).
  }{
    c = area(W)/sum[i] (1/lambda(x[i])).
  }
  This rescaling reduces the variability and bias of the estimate
  in small samples and in cases of very strong inhomogeneity.
  The default value of \code{normpower} is 1
  but the most sensible value is 2, which would correspond to rescaling
  the \code{lambda} values so that
  \eqn{
    \sum (1/\lambda(x_i)) = \mbox{area}(W).
  }{
    sum[i] (1/lambda(x[i])) = area(W).
  }
}
\value{
  A function value table (object of class \code{"fv"}).
  Essentially a data frame containing the variables
  \item{r}{
    the vector of values of the argument \eqn{r} 
    at which the inhomogeneous pair correlation function
    \eqn{g_{\rm inhom}(r)}{ginhom(r)} has been  estimated
  }
  \item{theo}{vector of values equal to 1,
    the theoretical value of \eqn{g_{\rm inhom}(r)}{ginhom(r)}
    for the Poisson process
  }
  \item{trans}{vector of values of \eqn{g_{\rm inhom}(r)}{ginhom(r)}
    estimated by translation correction
  }
  \item{iso}{vector of values of \eqn{g_{\rm inhom}(r)}{ginhom(r)}
    estimated by Ripley isotropic correction
  }
  as required.
}
\seealso{
  \code{\link{pcf}}, 
  \code{\link{pcf.ppp}}, 
  \code{\link{Kinhom}}
}
\examples{
  data(residualspaper)
  X <- residualspaper$Fig4b
  plot(pcfinhom(X, stoyan=0.2, sigma=0.1))
  fit <- ppm(X, ~polynom(x,y,2))
  plot(pcfinhom(X, lambda=fit, normpower=2))
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{nonparametric}
