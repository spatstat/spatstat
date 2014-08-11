\name{rMosaicField}
\alias{rMosaicField}
\title{Mosaic Random Field}
\description{
  Generate a realisation of a random field
  which is piecewise constant on the tiles of a given tessellation.
}
\usage{
rMosaicField(X, 
    rgen = function(n) { sample(0:1, n, replace = TRUE)},
    ...,
    rgenargs=NULL)
}
\arguments{
  \item{X}{
    A tessellation (object of class \code{"tess"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link{as.mask}}
    determining the pixel resolution.
  }
  \item{rgen}{
    Function that generates random values for the
    tiles of the tessellation.
  }
  \item{rgenargs}{
    List containing extra arguments that should be passed
    to \code{rgen} (typically specifying parameters of the
    distribution of the values).
  }
}
\details{
  This function generates a realisation of a random field
  which is piecewise constant on the tiles of the given
  tessellation \code{X}. The values in each tile
  are independent and identically distributed.
}
\value{
  A pixel image (object of class \code{"im"}).
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{rpoislinetess}},
  \code{\link{rMosaicSet}}
}
\examples{
   X <- rpoislinetess(3)
   plot(rMosaicField(X, runif))
   plot(rMosaicField(X, runif, dimyx=256))
   plot(rMosaicField(X, rnorm, rgenargs=list(mean=10, sd=2)))

   plot(rMosaicField(dirichlet(runifpoint(30)), rnorm))
}
\keyword{spatial}
\keyword{datagen}
