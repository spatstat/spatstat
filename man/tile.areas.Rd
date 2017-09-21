\name{tile.areas}
\alias{tile.areas}
\title{Compute Areas of Tiles in a Tessellation}
\description{
  Computes the area of each tile in a tessellation.
}
\usage{
tile.areas(x)
}
\arguments{
  \item{x}{A tessellation (object of class \code{"tess"}).}
}
\details{
  A tessellation is a collection of disjoint spatial regions
  (called \emph{tiles}) that fit together to form a larger spatial
  region. See \code{\link{tess}}.

  This command computes the area of each of the tiles 
  that make up the tessellation \code{x}.
  The result is a numeric vector
  in the same order as the tiles would be listed by \code{tiles(x)}.
}
\value{
  A numeric vector.
}
\seealso{
  \code{\link{tess}},
  \code{\link{tiles}},
  \code{\link{tilenames}},
  \code{\link{tiles.empty}}
}
\examples{
  A <- tess(xgrid=0:2,ygrid=0:2)
  tile.areas(A)
  v <- as.im(function(x,y){factor(round(x^2 + y^2))}, W=owin())
  E <- tess(image=v)
  tile.areas(E)
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{manip}
