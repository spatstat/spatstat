\name{grow.boxx}
\alias{grow.boxx}
\alias{grow.box3}
\title{Add margins to box in any dimension}
\description{
  Adds a margin to a box of class boxx.
}
\usage{
 grow.boxx(W, left, right = left)
 grow.box3(W, left, right = left)
}
\arguments{
  \item{W}{
    A box (object of class \code{"boxx"} or \code{"box3"}).
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
  Another object of the same class \code{"boxx"} or \code{"box3"}
  representing the window after margins are added.
}
\seealso{
  \code{\link{grow.rectangle}},
  \code{\link{boxx}},
  \code{\link{box3}}
}
\examples{
  w <- boxx(c(0,10), c(0,10), c(0,10), c(0,10))
  # add a margin of size 1 on both sides in all four dimensions
  b12 <- grow.boxx(w, 1)

  # add margin of size 2 at left, and margin of size 3 at right,
  # in each dimension.
  v <- grow.boxx(w, 2, 3)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{manip}
