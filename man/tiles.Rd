\name{tiles}
\alias{tiles}
\title{Extract List of Tiles in a Tessellation}
\description{
  Extracts a list of the tiles that make up a tessellation.
}
\usage{
tiles(x)
}
\arguments{
  \item{x}{A tessellation (object of class \code{"tess"}).}
}
\details{
  A tessellation is a collection of disjoint spatial regions
  (called \emph{tiles}) that fit together to form a larger spatial
  region. See \code{\link{tess}}.

  The tiles that make up the tessellation \code{x}
  are returned in a list.
}
\value{
  A list of windows (objects of class \code{"owin"}).
}
\seealso{
  \code{\link{tess}},
  \code{\link{tilenames}},
  \code{\link{tile.areas}},
  \code{\link{tiles.empty}}
}
\examples{
  A <- tess(xgrid=0:2,ygrid=0:2)
  tiles(A)
  v <- as.im(function(x,y){factor(round(x^2 + y^2))}, W=owin())
  E <- tess(image=v)
  tiles(E)
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{manip}
