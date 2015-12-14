\name{grow.boxx}
\alias{grow.boxx}
\title{Add margins to box in any dimension}
\description{
  Adds a margin to a box of class boxx.
}
\usage{
 grow.boxx(W, left, right = left)
}
\arguments{
  \item{W}{
    A box (object of class \code{"boxx"}).
  }
  \item{left}{Width of margin to be added to left endpoint
    of box side in every dimension.
    A single nonnegative number, or a vector of same length
    as the dimension of the box to add different left margin in each dimension.
  }
  \item{right}{Width of margin to be added to right endpoint
    of box side in every dimension.
    A single nonnegative number, or a vector of same length
    as the dimension of the box to add different right margin in each dimension.
  }
}
\value{
  Another object of class \code{"boxx"} representing the
  window after margins are added.
}
\seealso{
  \code{\link{grow.rectangle}},
  \code{\link{boxx}}
}
\examples{
  w <- boxx(c(0,10), c(0,10), c(0,10))
  # add a margin of size 1 on both sides in all three dimensions
  box11 <- grow.boxx(w, 1)

  # add margin of size 2 to the three left end points
  # and margin of size 3 to the right end points.
  v <- grow.boxx(w, c(2,2,2), c(3,3,3))
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}
  
}
\keyword{spatial}
\keyword{manip}
