\name{ippm}
\alias{ippm}
\title{
  Optimise Irregular Trend Parameters in Point Process Model
}
\description{
  Experimental extension to \code{ppm}.
  Find optimal values of the irregular trend parameters in a
  point process model using Fisher scoring algorithm.
}
\usage{
ippm(..., iScore=NULL,
          start=list(),
          covfunargs=start,
          maxiter=20, tol=1e-4, progress=TRUE, stepfactor=1,
          dbug=FALSE)
}
\arguments{
  \item{\dots}{
    Arguments passed to \code{\link[spatstat]{ppm}} to fit the point process model.
  }
  \item{iScore}{
    A named list of \R functions that compute the partial derivatives
    of \code{logf} with respect to each irregular parameter.
    See Details.
  }
  \item{start}{
    Named list containing initial values of the
    irregular parameters over which to optimise.
  }
  \item{covfunargs}{
    Argument passed to \code{\link[spatstat]{ppm}}.
    A named list containing values for \emph{all} irregular parameters
    required by the covariates in the model.
    Must include all the parameters named in \code{start}.
  }
  \item{maxiter}{
    Integer. Maximum number of iterations of Fisher scoring algorithm.
  }
  \item{tol}{
    Numeric value or vector. The algorithm stops when the
    difference between two successive estimates of the irregular
    parameter is less than \code{tol}.
  }
  \item{progress}{
    Logical. Whether to print progress reports.
  }
  \item{stepfactor}{
    Numeric value between 0 and 1 indicating that the change in the
    parameter between successive iterations is only a specified
    fraction of the step computed by the Newton-Raphson algorithm.
  }
  \item{dbug}{
    Logical. Whether to print debugging output.
  }
}
\details{
  This function is an experimental extension to the
  point process model fitting command \code{\link[spatstat]{ppm}}.
  The extension allows the trend of the model to include irregular parameters,
  which will be maximised by a Fisher scoring method. 

  For the sake of explanation,
  consider a Poisson point process with intensity function
  \eqn{\lambda(u)}{lambda(u)} at location \eqn{u}. Assume that
  \deqn{
    \lambda(u) = \exp(\alpha + \beta Z(u)) \, f(u, \gamma)
  }{
    lambda(u) = \exp(alpha + beta * Z(u)) * f(u, gamma)
  }
  where \eqn{\alpha,\beta,\gamma}{alpha,beta,gamma} are
  parameters to be estimated, \eqn{Z(u)} is a spatial covariate
  function, and \eqn{f} is some known function.
  Then the parameters
  \eqn{\alpha,\beta}{alpha,beta} are called \emph{regular} because they
  appear in a loglinear form; the parameter 
  \eqn{\gamma}{gamma} is called \emph{irregular}.
  
  To fit this model using \code{ippm}, we specify the
  intensity using the \code{trend} formula
  in the same way as usual for \code{\link[spatstat]{ppm}}.
  The trend formula is a representation of the log intensity.
  In the above example the log intensity is
  \deqn{
    \log\lambda(u) = \alpha + \beta Z(u) + \log f(u, \gamma)
  }{
    log(lambda(u)) =  alpha + beta * Z(u) + log(f(u, gamma))
  }
  So the model above would be encoded with the trend formula
  \code{~Z + offset(log(f))}. Note that the irregular part of the model
  is an \emph{offset} term, which means that it is included in the log trend
  as it is, without being multiplied by another regular parameter.

  To perform Fisher scoring we also need the derivative
  of \eqn{\log f(u,\gamma)}{log(f(u,gamma))} with
  respect to \eqn{\gamma}{gamma}. We call this the
  \emph{irregular score}. The user must write an \R function
  that computes the irregular score for any value of
  \eqn{\gamma}{gamma} at any location \code{(x,y)}.
  
  Thus, to code such a problem,
  \enumerate{
    \item The argument \code{trend} should define the
    log intensity, with the irregular part as an offset;
    \item The argument \code{start} should be a list
    containing initial values of each of the irregular parameters;
    \item The argument \code{iScore} must be a list (with one entry
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
    f((x,y), (gamma, delta)) = 1 + \exp(gamma - delta * x^3)
  }

  Arguments \code{\dots} passed to \code{\link[spatstat]{ppm}} may
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
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link[spatstat]{ppm}}
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
  ippm(X, ~Z + offset(log(f)),
       covariates=list(Z=Z, f=f),
       iScore=Dlogf,
       start=list(gamma=1, delta=1),
       tol=0.01, nd=nd)
}
\keyword{spatial}
\keyword{models}
