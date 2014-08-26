\name{domain}
\alias{domain}
\alias{domain.ppp}
\alias{domain.psp}
\alias{domain.im}
\alias{domain.ppx}
\alias{domain.pp3}
\alias{domain.lpp}
\alias{domain.ppm}
\alias{domain.kppm}
\alias{domain.lpp}
\alias{domain.lppm}
\alias{domain.msr}
\alias{domain.quad}
\alias{domain.quadratcount}
\alias{domain.quadrattest}
\alias{domain.tess}
\alias{domain.im}
\alias{domain.layered}
\alias{domain.distfun}
\alias{domain.nnfun}
\alias{domain.funxy}
\alias{domain.rmhmodel}
\alias{domain.leverage.ppm}
\alias{domain.influence.ppm}
\title{
  Extract the Domain of any Spatial Object
}
\description{
  Given a spatial object such as a point pattern, in any number of dimensions,
  this function extracts the spatial domain in which the object is defined.
}
\usage{
   domain(X, \dots)

   \method{domain}{ppp}(X, \dots)

   \method{domain}{psp}(X, \dots)

   \method{domain}{im}(X, \dots)

   \method{domain}{ppx}(X, \dots)

   \method{domain}{pp3}(X, \dots)

   \method{domain}{lpp}(X, \dots)

 \method{domain}{ppm}(X, \dots, from=c("points", "covariates"))

 \method{domain}{kppm}(X, \dots, from=c("points", "covariates"))

 \method{domain}{lpp}(X, \dots)

 \method{domain}{lppm}(X, \dots)

 \method{domain}{msr}(X, \dots)

 \method{domain}{quad}(X, \dots)

 \method{domain}{quadratcount}(X, \dots)

 \method{domain}{quadrattest}(X, \dots)

 \method{domain}{tess}(X, \dots)

 \method{domain}{im}(X, \dots)

 \method{domain}{layered}(X, \dots)

 \method{domain}{distfun}(X, \dots)

 \method{domain}{nnfun}(X, \dots)

 \method{domain}{funxy}(X, \dots)

 \method{domain}{rmhmodel}(X, \dots)

 \method{domain}{leverage.ppm}(X, \dots)

 \method{domain}{influence.ppm}(X, \dots)
}
\arguments{
  \item{X}{
    A spatial object such as a point pattern (in any number
    of dimensions), line segment pattern or pixel image.
  }
  \item{\dots}{
    Extra arguments. They are ignored by all the methods listed here.
  }
  \item{from}{Character string. See Details.}
}
\details{
  The function \code{domain} is generic.

  For a spatial object \code{X} in any number of dimensions, 
  \code{domain(X)} extracts the spatial domain in which \code{X} is
  defined.

  For a two-dimensional object \code{X}, typically \code{domain(X)}
  is the same as \code{domain(X)}.

  The exception is that, if \code{X} is a point pattern on a linear network
  (class \code{"lpp"}) or a point process model on a linear network
  (class \code{"lppm"}), then \code{domain(X)} is the linear network
  on which the points lie, while \code{Window(X)} is the two-dimensional
  window containing the linear network.

  The argument \code{from} applies when \code{X} is a fitted 
  point process model
  (object of class \code{"ppm"} or \code{"kppm"}).
  If \code{from="data"} (the default),
  \code{domain} extracts the window of the original point
  pattern data to which the model was fitted.
  If \code{from="covariates"} then \code{domain} returns the
  window in which the spatial covariates of the model were provided.
}
\value{
  A spatial object representing the domain of \code{X}.
  Typically a window (object of class \code{"owin"}),
  a three-dimensional box (\code{"box3"}), a multidimensional
  box (\code{"boxx"}) or a linear network (\code{"linnet"}).
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\seealso{
  \code{\link{domain}}
}
\examples{
  domain(cells)
  domain(bei.extra$elev)
  domain(chicago)
}
\keyword{spatial}
\keyword{manip}
