\name{methods.rho2hat}
\alias{methods.rho2hat} %DoNotExport
\alias{print.rho2hat}
\alias{plot.rho2hat}
\title{
  Methods for Intensity Functions of Two Spatial Covariates
}
\description{
  These are methods for the class \code{"rho2hat"}. 
}
\usage{
\method{plot}{rho2hat}(x, \dots, do.points=FALSE)
\method{print}{rho2hat}(x, \dots)
}
\arguments{
  \item{x}{
    An object of class \code{"rho2hat"}.
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
  \item{do.points}{
    Logical value indicating whether to plot the observed values of
    the covariates at the data points.
  }
}
\details{
  These functions are methods for the generic commands
  \code{\link{print}} and
  \code{\link{plot}} 
  for the class \code{"rho2hat"}.

  An object of class \code{"rho2hat"} is an estimate
  of the intensity of a point process, as a function of two
  given spatial covariates. See \code{\link{rho2hat}}.
  
  The method \code{plot.rho2hat} displays the estimated function
  \eqn{\rho}{rho} using \code{\link{plot.fv}}, and optionally
  adds a \code{\link{rug}} plot of the observed values of the covariate.
}
\value{
  \code{NULL}.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\seealso{
  \code{\link{rho2hat}}
}
\examples{
  data(bei)
  attach(bei.extra)
  r2 <- rho2hat(bei, elev, grad)
  r2
  plot(r2) 
}
\keyword{spatial}
\keyword{methods}
