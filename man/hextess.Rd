\name{hextess}
\alias{hexgrid}
\alias{hextess}
\title{
  Hexagonal Grid or Tessellation
}
\description{
  Construct a hexagonal grid of points,
  or a hexagonal tessellation.
}
\usage{
hexgrid(W, s, offset = c(0, 0), origin=NULL, trim = TRUE)

hextess(W, s, offset = c(0, 0), origin=NULL, trim = TRUE)
}

\arguments{
  \item{W}{
    Window in which to construct the hexagonal grid or tessellation.
    An object of class \code{"owin"}.
  }
  \item{s}{
    Side length of hexagons. A positive number.
  }
  \item{offset}{
    Numeric vector of length 2 specifying a shift of the
    hexagonal grid. See Details.
  }
  \item{origin}{
    Numeric vector of length 2 specifying the initial origin
    of the hexagonal grid, before the offset is applied.
    See Details.
  }
  \item{trim}{
    Logical value indicating whether to restrict the result to
    the window \code{W}. See Details.
  }
}
\details{
  \code{hexgrid} constructs a hexagonal grid of points
  on the window \code{W}. If \code{trim=TRUE} (the default),
  the grid is intersected with \code{W} so that all points lie
  inside \code{W}. If \code{trim=FALSE}, then we retain all grid points
  which are the centres of hexagons that intersect \code{W}.

  \code{hextess} constructs a tessellation of hexagons
  on the window \code{W}. If \code{trim=TRUE} (the default),
  the tessellation is restricted to the interior of \code{W},
  so that there will be some fragmentary hexagons near the
  boundary of \code{W}. If \code{trim=FALSE}, the tessellation
  consists of all hexagons which intersect \code{W}.

  The points of \code{hexgrid(...)} are the
  centres of the tiles of \code{hextess(...)} 
  in the same order.
  
  In the initial position of the grid or tessellation,
  one of the grid points (tile centres) is placed at the 
  \code{origin}, which defaults to the midpoint of the
  bounding rectangle of \code{W}. The grid can be shifted
  relative to this origin by specifing the \code{offset}.
}
\value{
  The value of \code{hexgrid} is a point pattern (object of class
  \code{"ppp"}).

  The value of \code{hextess} is a tessellation (object of class \code{"tess"}).
}
\seealso{
  \code{\link{tess}}
}
\examples{
  if(interactive()) {
    W <- Window(chorley)
    s <- 0.7
  } else {
    W <- letterR
    s <- 0.3
  }
  plot(hextess(W, s))
  plot(hexgrid(W, s), add=TRUE)
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{datagen}
