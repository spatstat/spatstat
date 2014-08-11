\name{rMosaicSet}
\alias{rMosaicSet}
\title{Mosaic Random Set}
\description{
  Generate a random set by taking a random selection of 
  tiles of a given tessellation.  
}
\usage{
rMosaicSet(X, p=0.5)
}
\arguments{
  \item{X}{
    A tessellation (object of class \code{"tess"}).
  }
  \item{p}{
    Probability of including a given tile.
    A number strictly between 0 and 1.
  }
}
\details{
  Given a tessellation \code{X}, this function 
  randomly selects some of the tiles of \code{X},
  including each tile with probability \eqn{p} independently of
  the other tiles. The selected tiles are then combined to form a
  set in the plane. 

  One application of this is Switzer's (1965) example of a random set
  which has a Markov property. It is constructed by generating \code{X}
  according to a Poisson line tessellation (see \code{\link{rpoislinetess}}).
}
\value{
  A window (object of class \code{"owin"}).
}
\references{
  Switzer, P.
  A random set process in the plane with a Markovian property.
  \emph{Annals of Mathematical Statistics} \bold{36} (1965) 1859--1863.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{rpoislinetess}},
  \code{\link{rMosaicField}}
}
\examples{
   # Switzer's random set
   X <- rpoislinetess(3)
   plot(rMosaicSet(X, 0.5), col="green", border=NA)

   # another example
   plot(rMosaicSet(dirichlet(runifpoint(30)), 0.4))
}
\keyword{spatial}
\keyword{datagen}
