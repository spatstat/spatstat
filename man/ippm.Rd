\name{ippm}
\alias{ippm}
\title{
  Fit Point Process Model Involving Irregular Trend Parameters
}
\description{
  Experimental extension to \code{ppm}
  which finds optimal values of the irregular trend parameters in a
  point process model.
}
\usage{
ippm(Q, \dots,
          iScore=NULL,
          start=list(),
          covfunargs=start,
          nlm.args=list(stepmax=1/2),
          silent=FALSE,
          warn.unused=TRUE)
}
\arguments{
  \item{Q,\dots}{
    Arguments passed to \code{\link{ppm}}
    to fit the point process model.
  }
  \item{iScore}{
    Optional. A named list of \R functions that compute the partial derivatives
    of the logarithm of the trend, with respect to each irregular parameter.
    See Details.
  }
  \item{start}{
    Named list containing initial values of the
    irregular parameters over which to optimise.
  }
  \item{covfunargs}{
    Argument passed to \code{\link{ppm}}.
    A named list containing values for \emph{all} irregular parameters
    required by the covariates in the model.
    Must include all the parameters named in \code{start}.
  }
  \item{nlm.args}{
    Optional list of arguments passed to \code{\link[stats]{nlm}}
    to control the optimization algorithm.
  }
  \item{silent}{
    Logical. Whether to print warnings if the optimization algorithm
    fails to converge.
  }
  \item{warn.unused}{
    Logical. Whether to print a warning if some of the parameters
    in \code{start} are not used in the model.
  }
}
\details{
  This function is an experimental extension to the
  point process model fitting command \code{\link{ppm}}.
  The extension allows the trend of the model to include irregular parameters,
  which will be maximised by a Newton-type iterative
  method, using \code{\link[stats]{nlm}}.

  For the sake of explanation,
  consider a Poisson point process with intensity function
  \eqn{\lambda(u)}{\lambda(u)} at location \eqn{u}. Assume that
  \deqn{
    \lambda(u) = \exp(\alpha + \beta Z(u)) \, f(u, \gamma)
  }{
    \lambda(u) = exp(\alpha + \beta * Z(u)) * f(u, \gamma)
  }
  where \eqn{\alpha,\beta,\gamma} are
  parameters to be estimated, \eqn{Z(u)} is a spatial covariate
  function, and \eqn{f} is some known function.
  Then the parameters
  \eqn{\alpha,\beta} are called \emph{regular} because they
  appear in a loglinear form; the parameter 
  \eqn{\gamma} is called \emph{irregular}.
  
  To fit this model using \code{ippm}, we specify the
  intensity using the \code{trend} formula
  in the same way as usual for \code{\link{ppm}}.
  The trend formula is a representation of the log intensity.
  In the above example the log intensity is
  \deqn{
    \log\lambda(u) = \alpha + \beta Z(u) + \log f(u, \gamma)
  }{
    log(\lambda(u)) =  \alpha + \beta * Z(u) + log(f(u, \gamma))
  }
  So the model above would be encoded with the trend formula
  \code{~Z + offset(log(f))}. Note that the irregular part of the model
  is an \emph{offset} term, which means that it is included in the log trend
  as it is, without being multiplied by another regular parameter.

  The optimisation runs faster if we specify the derivative
  of \eqn{\log f(u,\gamma)}{log(f(u,\gamma))} with
  respect to \eqn{\gamma}. We call this the
  \emph{irregular score}. To specify this, the user must write an \R function
  that computes the irregular score for any value of
  \eqn{\gamma} at any location \code{(x,y)}.
  
  Thus, to code such a problem,
  \enumerate{
    \item The argument \code{trend} should define the
    log intensity, with the irregular part as an offset;
    \item The argument \code{start} should be a list
    containing initial values of each of the irregular parameters;
    \item The argument \code{iScore}, if provided,
    must be a list (with one entry
    for each entry of \code{start}) of functions
    with arguments \code{x,y,\dots}, that evaluate the partial derivatives
    of \eqn{\log f(u,\gamma)}{log(f(u,gamma))} with
    respect to each irregular parameter.
  }
  
  The coded example below illustrates the model with two irregular
  parameters \eqn{\gamma,\delta}{gamma,delta} and irregular term
  \deqn{
    f((x,y), (\gamma, \delta)) = 1 + \exp(\gamma - \delta x^3)
  }{
    f((x,y), (\gamma, \delta)) = 1 + \exp(\gamma - \delta * x^3)
  }

  Arguments \code{\dots} passed to \code{\link{ppm}} may
  also include \code{interaction}. In this case the model is not
  a Poisson point process but a more general Gibbs point process;
  the trend formula \code{trend} 
  determines the first-order trend
  of the model (the first order component of the conditional intensity),
  not the intensity.
}
\value{
  A fitted point process model (object of class \code{"ppm"}).
}
\author{\spatstatAuthors.}
\seealso{
  \code{\link{ppm}},
  \code{\link{profilepl}}
}
\examples{
  nd <- 32
  \testonly{nd <- 10}
  
  gamma0 <- 3
  delta0 <- 5
  POW <- 3
  # Terms in intensity
  Z <- function(x,y) { -2*y }
  f <- function(x,y,gamma,delta) { 1 + exp(gamma - delta * x^POW) }
  # True intensity
  lamb <- function(x,y,gamma,delta) { 200 * exp(Z(x,y)) * f(x,y,gamma,delta) }
  # Simulate realisation
  lmax <- max(lamb(0,0,gamma0,delta0), lamb(1,1,gamma0,delta0))
  set.seed(42)
  X <- rpoispp(lamb, lmax=lmax, win=owin(), gamma=gamma0, delta=delta0)
  # Partial derivatives of log f
  DlogfDgamma <- function(x,y, gamma, delta) {
    topbit <- exp(gamma - delta * x^POW)
    topbit/(1 + topbit)
  }
  DlogfDdelta <- function(x,y, gamma, delta) {
    topbit <- exp(gamma - delta * x^POW)
    - (x^POW) * topbit/(1 + topbit)
  }
  # irregular score
  Dlogf <- list(gamma=DlogfDgamma, delta=DlogfDdelta)
  # fit model
  ippm(X ~Z + offset(log(f)),
       covariates=list(Z=Z, f=f),
       iScore=Dlogf,
       start=list(gamma=1, delta=1),
       nd=nd)
}
\keyword{spatial}
\keyword{models}
