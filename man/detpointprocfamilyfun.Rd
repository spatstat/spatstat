\name{detpointprocfamilyfun}
\alias{detpointprocfamilyfun}
\title{Construct a New Determinantal Point Process Model Family Function}
\description{
  Function to ease the implementation of a new determinantal
  point process model family.
}
\usage{detpointprocfamilyfun(kernel = NULL,
    specden = NULL, basis = "fourierbasis", 
    convkernel = NULL, Kfun = NULL, valid = NULL, intensity = NULL, 
    dim = 2, name = "User-defined", isotropic = TRUE, range = NULL, 
    parbounds = NULL, specdenrange = NULL, startpar = NULL, \dots)
}
\arguments{
  \item{kernel}{
    function specifying the kernel.
    May be set to \code{NULL}. See Details.
  }
  \item{specden}{
    function specifying the spectral density.
    May be set to \code{NULL}. See Details.
  }
  \item{basis}{
    character string giving the name of the basis.
    Defaults to the Fourier basis. See Details.
  }
  \item{convkernel}{
    function specifying the k-fold auto-convolution of the kernel.
    May be set to \code{NULL}. See Details.
  }
  \item{Kfun}{
    function specifying the K-function.
    May be set to \code{NULL}. See Details.
  }
  \item{valid}{
    function determining whether a given set of parameter values
    yields a valid model. May be set to \code{NULL}. See Examples.
  }
  \item{intensity}{
    character string specifying which parameter is the intensity in the
    model family. Should be \code{NULL} if the model family has no intensity
    parameter.
  }
  \item{dim}{
    character strig specifying which parameter is the dimension of the
    state space in this model family (if any).
    Alternatively a positive integer specifying the dimension.
  }
  \item{name}{
    character string giving the name of the model family used for
    printing.
  }
  \item{isotropic}{
    logical value indicating whether or not the model is isotropic.
  }
  \item{range}{
    function determining the interaction range of the model. May be
    set to \code{NULL}. See Examples.
  }
  \item{parbounds}{
    function determining the bounds for each model parameter when all
    other parameters are fixed. May be set to \code{NULL}. See Examples.
  }
  \item{specdenrange}{
    function specifying the the range of the spectral density if it
    is finite (only the case for very few models). May be set to
    \code{NULL}.
  }
  \item{startpar}{
    function determining starting values for parameters in any estimation
    algorithm. May be set to \code{NULL}. See Examples.
  }
  \item{\dots}{
    Additional arguments for inclusion in the returned model object. These
    are not checked in any way.
  }
}
\details{
  A determinantal point process family is specified either
  in terms of a kernel (a positive semi-definite function, i.e. a
  covariance function) or a spectral density, or preferably both. 
  One of these can be \code{NULL} if it is unknown, but not both.
  When both are supplied they must have the same arguments.
  The first argument gives
  the values at which the function should be evaluated. In general
  the function should accept an  \eqn{n} by \eqn{d} matrix or
  \code{data.frame} specifying \eqn{n (>=0)}
  points in dimension \eqn{d}. If the model is isotropic it only needs to
  accept a non-negative valued numeric of length \eqn{n}. (In fact
  there is currently almost no support for non-isotropic models, so
  it is recommended not to specify such a model.) The name
  of this argument could be chosen freely, but \eqn{x} is
  recommended. The remaining arguments are the parameters of the
  model. If one of these is an intensity parameter the name should
  be mentioned in the argument \code{intensity}. If one of these
  specifies the dimension of the model it should be mentioned in the
  argument \code{dim}.

  The kernel and spectral density is with respect to a
  specific set of basis functions, which is typically the Fourier
  basis. However this can be changed to any user-supplied basis in
  the argument \code{basis}. If such an alternative is supplied it
  must be the name of a function expecting the same arguments as
  \code{\link{fourierbasis}} and returning the results in the same
  form as \code{\link{fourierbasis}}.

  If supplied, the arguments of convkernel must obey the
  following: first argument should be like the first argument of
  kernel and/or specden (see above). The second argument (preferably
  called \code{k}) should be the positive integer specifying how many times
  the auto-convolution is done (i.e. the \eqn{k} in \eqn{k}-fold
  auto-convolution). The remaining arguments must agree with the
  arguments of \code{kernel} and/or \code{specden} (see above).

  If supplied, the arguments of \code{Kfun} should be like the
  arguments of \code{kernel} and \code{specden} (see above).
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\examples{
  ## Example of how to define the Gauss family
exGauss <- detpointprocfamilyfun(
    name="Gaussian",
    kernel=function(x, lambda, alpha, d){
        lambda*exp(-(x/alpha)^2)
    },
    specden=function(x, lambda, alpha, d){
        lambda * (sqrt(pi)*alpha)^d * exp(-(x*alpha*pi)^2)
    },
    convkernel=function(x, k, lambda, alpha, d){
        logres <- k*log(lambda*pi*alpha^2) - log(pi*k*alpha^2) - x^2/(k*alpha^2)
        return(exp(logres))
    },
    Kfun = function(x, lambda, alpha, d){
        pi*x^2 - pi*alpha^2/2*(1-exp(-2*x^2/alpha^2))
    },
    valid=function(lambda, alpha, d){
        lambda>0 && alpha>0 && d>=1 && lambda <= (sqrt(pi)*alpha)^(-d)
    },
    isotropic=TRUE,
    intensity="lambda",
    dim="d",
    range=function(alpha, bound = .99){
        if(missing(alpha))
            stop("The parameter alpha is missing.")
        if(!(is.numeric(bound)&&bound>0&&bound<1))
            stop("Argument bound must be a numeric between 0 and 1.")
        return(alpha*sqrt(-log(sqrt(1-bound))))
    },
    parbounds=function(name, lambda, alpha, d){
        switch(name,
               lambda = c(0, (sqrt(pi)*alpha)^(-d)),
               alpha = c(0, lambda^(-1/d)/sqrt(pi)),
               stop("Parameter name misspecified")
               )
    },
    startpar=function(model, X){
        rslt <- NULL
        if("lambda" \%in\% model$freepar){
            lambda <- intensity(X)
            rslt <- c(rslt, "lambda" = lambda)
            model <- update(model, lambda=lambda)
        }
        if("alpha" \%in\% model$freepar){
            alpha <- .8*dppparbounds(model, "alpha")[2]
            rslt <- c(rslt, "alpha" = alpha)
        }
        return(rslt)
    }
    )
  exGauss
  m <- exGauss(lambda=100, alpha=.05, d=2)
  m
}
\keyword{spatial}
\keyword{models}
