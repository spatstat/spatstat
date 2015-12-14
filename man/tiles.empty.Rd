\name{tiles.empty}
\alias{tiles.empty}
\title{Check For Empty Tiles in a Tessellation}
\description{
  Checks whether each tile in a tessellation is empty or non-empty.
}
\usage{
tiles.empty(x)
}
\arguments{
  \item{x}{A tessellation (object of class \code{"tess"}).}
}
\details{
  A tessellation is a collection of disjoint spatial regions
  (called \emph{tiles}) that fit together to form a larger spatial
  region. See \code{\link{tess}}.

  It is possible for some tiles of a tessellation to be empty.
  For example, this can happen
  when the tessellation \code{x} is obtained by restricting
  another tessellation \code{y} to a smaller spatial domain \code{w}.

  The function \code{tiles.empty} checks whether each tile is empty
  or non-empty. The result is a logical vector,
  with entries equal to \code{TRUE} when the corresponding tile is
  empty. Results are given in the same order
  as the tiles would be listed by \code{tiles(x)}.
}
\value{
  A logical vector.
}
\seealso{
  \code{\link{tess}},
  \code{\link{tiles}},
  \code{\link{tilenames}},
  \code{\link{tile.areas}}
}
\examples{
  A <- tess(xgrid=0:2,ygrid=0:2)
  tiles.empty(A)
  v <- as.im(function(x,y){factor(round(x^2 + y^2))}, W=owin())
  E <- tess(image=v)
  tiles.empty(E)
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}
  
}
\keyword{spatial}
\keyword{manip}
