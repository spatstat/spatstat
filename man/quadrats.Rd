\name{quadrats}
\alias{quadrats}
\title{Divide Region into Quadrats}
\description{
  Divides window into rectangular quadrats and 
  returns the quadrats as a tessellation.
}
\usage{
quadrats(X, nx = 5, ny = nx, xbreaks = NULL, ybreaks = NULL, keepempty=FALSE)
}
\arguments{
  \item{X}{
    A window (object of class \code{"owin"})
    or anything that can be coerced to a window using
    \code{\link{as.owin}}, such as a point pattern.
  }
  \item{nx,ny}{
    Numbers of quadrats in the \eqn{x} and \eqn{y} directions.
    Incompatible with \code{xbreaks} and \code{ybreaks}.
  }
  \item{xbreaks}{
    Numeric vector giving the \eqn{x} coordinates of the
    boundaries of the quadrats. Incompatible with \code{nx}.
  }
  \item{ybreaks}{
    Numeric vector giving the \eqn{y} coordinates of the
    boundaries of the quadrats. Incompatible with \code{ny}.
  }
  \item{keepempty}{
    Logical value indicating whether to delete or retain
    empty quadrats. See Details.
  }
}
\details{
  If the window \code{X} is a rectangle, it is divided into
  an \code{nx * ny} grid of rectangular tiles or `quadrats'.

  If \code{X} is not a rectangle, then the bounding rectangle of
  \code{X} is first divided into an \code{nx * ny} grid of rectangular
  tiles, and these tiles are then intersected with the window \code{X}.

  The resulting tiles are returned as a tessellation (object of class
  \code{"tess"}) which can be plotted and used in other analyses.

  If \code{xbreaks} is given, it should be a numeric vector
  giving the \eqn{x} coordinates of the quadrat boundaries.
  If it is not given, it defaults to a
  sequence of \code{nx+1} values equally spaced
  over the range of \eqn{x} coordinates in the window \code{Window(X)}.

  Similarly if \code{ybreaks} is given, it should be a numeric
  vector giving the \eqn{y} coordinates of the quadrat boundaries.
  It defaults to a vector of \code{ny+1} values
  equally spaced over the range of \eqn{y} coordinates in the window.
  The lengths of \code{xbreaks} and \code{ybreaks} may be different.

  By default (if \code{keepempty=FALSE}), any rectangular tile which
  does not intersect the window \code{X} is
  ignored, and only the non-empty intersections are treated as quadrats,
  so the tessellation may consist of fewer than \code{nx * ny} tiles.
  If \code{keepempty=TRUE}, empty intersections are retained,
  and the tessellation always contains exactly \code{nx * ny} tiles,
  some of which may be empty.
}
\value{
  A tessellation (object of class \code{"tess"}) as described under
  \code{\link{tess}}.
}
\examples{
 W <- square(10)
 Z <- quadrats(W, 4, 5)
 plot(Z)

 data(letterR)
 plot(quadrats(letterR, 5, 7))
}
\seealso{
  For calculations using quadrats, see
  \code{\link{quadratcount}},
  \code{\link{quadrat.test}},
  \code{\link{quadratresample}}

  For other kinds of tessellations, see 
  \code{\link{tess}},
  \code{\link{hextess}},
  \code{\link{venn.tess}},
  \code{\link{polartess}},
  \code{\link{dirichlet}}, \code{\link{delaunay}},
  \code{\link{rpoislinetess}} and \code{\link{quantess}}.

}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{utilities}
\keyword{datagen}
