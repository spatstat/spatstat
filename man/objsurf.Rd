\name{objsurf}
\alias{objsurf}
\alias{objsurf.kppm}
\alias{objsurf.minconfit}
\title{
  Objective Function Surface
}
\description{
  For a model that was fitted by optimisation,
  compute the values of the objective function in a neighbourhood
  of the optimal value.
}
\usage{
objsurf(x, \dots)

\method{objsurf}{kppm}(x, ..., ngrid = 32, ratio = 1.5, verbose = TRUE)

\method{objsurf}{minconfit}(x, ..., ngrid = 32, ratio = 1.5, verbose = TRUE)
}
\arguments{
  \item{x}{
    Some kind of model that was fitted
    by finding the optimal value of an objective function. 
    An object of class \code{"kppm"} or \code{"minconfit"}.
  }
  \item{\dots}{
    Extra arguments are usually ignored.
  }
  \item{ngrid}{
    Number of grid points to evaluate along each axis.
    Either a single integer, or a pair of integers.
    For example \code{ngrid=32} would mean a \code{32 * 32} grid.
  }
  \item{ratio}{
    Number greater than 1 determining the range of parameter values
    to be considered. If the optimal parameter value is \code{opt}
    then the objective function will be evaluated for
    values between \code{opt/ratio} and \code{opt * ratio}.
  }
  \item{verbose}{
    Logical value indicating whether to print progress reports.
  }
}
\details{
  The object \code{x} should be some kind of model that was fitted
  by maximising or minimising the value of an objective function.
  The objective function will be evaluated on a grid of
  values of the model parameters. 

  Currently the following types of objects are accepted:
  \itemize{
    \item
    an object of class \code{"kppm"} representing a
    cluster point process or Cox point process. 
    See \code{\link{kppm}}.
    \item
    an object of class \code{"minconfit"} representing a
    minimum-contrast fit between a summary function and its
    theoretical counterpart. 
    See \code{\link{mincontrast}}.
  }
  The result is an object of class \code{"objsurf"} which can be
  printed and plotted: see \code{\link{methods.objsurf}}.
}
\value{
  An object of class \code{"objsurf"} which can be
  printed and plotted.
  Essentially a list containing entries \code{x}, \code{y}, \code{z}
  giving the parameter values and objective function values.
}
\author{
  Adrian Baddeley and Ege Rubak.
}
\seealso{
  \code{\link{methods.objsurf}},
  \code{\link{kppm}}, 
  \code{\link{mincontrast}}
}
\examples{
   fit <- kppm(redwood, ~1, "Thomas")
   os <- objsurf(fit)

   if(interactive()) {
     plot(os)
     contour(os, add=TRUE)
     persp(os)
   }
}
\keyword{spatial}
\keyword{models}
