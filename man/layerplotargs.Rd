\name{layerplotargs}
\alias{layerplotargs}
\alias{layerplotargs<-}
\title{
  Extract or Replace the Plot Arguments of a Layered Object
}
\description{
  Extracts or replaces the plot arguments of a layered object.
}
\usage{
layerplotargs(L)

layerplotargs(L) <- value
}
\arguments{
  \item{L}{
    An object of class \code{"layered"}
    created by the function \code{\link{layered}}.
  }
  \item{value}{
    Replacement value. 
    A list, with the same length as \code{L},
    whose elements are lists of plot arguments.
  }
}
\details{
  These commands extract or replace the \code{plotargs}
  in a layered object. See \code{\link{layered}}.
}
\value{
  \code{layerplotargs} returns a list of lists of plot arguments.

  \code{"layerplotargs<-"} returns the updated object 
  of class \code{"layered"}.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{layered}},
  \code{\link{methods.layered}},
  \code{\link{[.layered}}.
}
\examples{
   W <- square(2)
   L <- layered(W=W, X=cells)
   layerplotargs(L)$X <- list(pch=16)
}
\keyword{spatial}
\keyword{hplot}
