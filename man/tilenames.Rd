\name{tilenames}
\alias{tilenames}
\alias{tilenames<-}
\alias{tilenames.tess}
\alias{tilenames<-.tess}
\alias{tilenames.lintess}
\alias{tilenames<-.lintess}
\title{Names of Tiles in a Tessellation}
\description{
  Extract or Change the Names of the Tiles in a Tessellation.
}
\usage{
tilenames(x)
tilenames(x) <- value

\method{tilenames}{tess}(x)

\method{tilenames}{lintess}(x)

\method{tilenames}{tess}(x) <- value

\method{tilenames}{lintess}(x) <- value
}
\arguments{
  \item{x}{
    A tessellation (object of class \code{"tess"})
    or a tessellation on a linear network (object of class \code{"lintess"}).
  }
  \item{value}{Character vector giving new names for the tiles.}
}
\details{
  These functions extract or change the names of the 
  tiles that make up the tessellation \code{x}.

  If the tessellation is a regular grid, the tile names
  cannot be changed.
}
\value{
  \code{tilenames} returns a character vector.
}
\seealso{
  \code{\link{tess}}, \code{\link{lintess}},
  \code{\link{tiles}}
}
\examples{
  D <- dirichlet(runifpoint(10))
  tilenames(D)
  tilenames(D) <- paste("Cell", 1:10)
  tilenames(D)

  B <- lineardirichlet(runiflpp(5, simplenet))
  tilenames(B)
  tilenames(B) <- letters[1:5]
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{manip}
