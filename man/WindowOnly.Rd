\name{WindowOnly}
\alias{Window.ppm}
\alias{Window.kppm}
\alias{Window.dppm}
\alias{Window.lpp}
\alias{Window.lppm}
\alias{Window.msr}
\alias{Window.quad}
\alias{Window.quadratcount}
\alias{Window.quadrattest}
\alias{Window.tess}
\alias{Window.layered}
\alias{Window.distfun}
\alias{Window.nnfun}
\alias{Window.funxy}
\alias{Window.rmhmodel}
\alias{Window.leverage.ppm}
\alias{Window.influence.ppm}
\title{Extract Window of Spatial Object}
\description{
  Given a spatial object (such as a point pattern or pixel image)
  in two dimensions, these functions extract the window
  in which the object is defined.
}
\usage{
 \method{Window}{ppm}(X, \dots, from=c("points", "covariates"))

 \method{Window}{kppm}(X, \dots, from=c("points", "covariates"))

 \method{Window}{dppm}(X, \dots, from=c("points", "covariates"))

 \method{Window}{lpp}(X, \dots)

 \method{Window}{lppm}(X, \dots)

 \method{Window}{msr}(X, \dots)

 \method{Window}{quad}(X, \dots)

 \method{Window}{quadratcount}(X, \dots)

 \method{Window}{quadrattest}(X, \dots)

 \method{Window}{tess}(X, \dots)

 \method{Window}{layered}(X, \dots)

 \method{Window}{distfun}(X, \dots)

 \method{Window}{nnfun}(X, \dots)

 \method{Window}{funxy}(X, \dots)

 \method{Window}{rmhmodel}(X, \dots)

 \method{Window}{leverage.ppm}(X, \dots)

 \method{Window}{influence.ppm}(X, \dots)
}
\arguments{
  \item{X}{A spatial object.}
  \item{\dots}{Ignored.}
  \item{from}{Character string. See Details.}
}
\value{
  An object of class \code{"owin"} (see \code{\link{owin.object}})
  specifying an observation window.
}
\details{
  These are methods for the generic function \code{\link{Window}}
  which extract the spatial window in which the object \code{X}
  is defined. 

  The argument \code{from} applies when \code{X} is a fitted 
  point process model
  (object of class \code{"ppm"}, \code{"kppm"} or \code{"dppm"}).
  If \code{from="data"} (the default),
  \code{Window} extracts the window of the original point
  pattern data to which the model was fitted.
  If \code{from="covariates"} then \code{Window} returns the
  window in which the spatial covariates of the model were provided.
}
\seealso{
  \code{\link{Window}},
  \code{\link{Window.ppp}},
  \code{\link{Window.psp}}.
  
  \code{\link{owin.object}}
}
\examples{
   X <- quadratcount(cells, 4)
   Window(X)
}
\author{\adrian
  
  
  \rolf
  
  and \ege 
}
\keyword{spatial}
\keyword{manip}
 
 
