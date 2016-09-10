\name{grow.rectangle}
\alias{grow.rectangle}
\title{Add margins to rectangle}
\description{
  Adds a margin to a rectangle.
}
\usage{
 grow.rectangle(W, xmargin=0, ymargin=xmargin, fraction=NULL)
}
\arguments{
  \item{W}{
    A window (object of class \code{"owin"}).
    Must be of type \code{"rectangle"}.
  }
  \item{xmargin}{Width of horizontal margin to be added.
    A single nonnegative number, or a vector of length 2
    indicating margins of unequal width at left and right.
  }
  \item{ymargin}{Height of vertical margin to be added.
    A single nonnegative number, or a vector of length 2
    indicating margins of unequal width at bottom and top.
  }
  \item{fraction}{
    Fraction of width and height to be added.
    A number greater than zero, or a numeric vector of length 2 indicating
    different fractions of width and of height, respectively.
    Incompatible with specifying \code{xmargin} and \code{ymargin}.
  }
}
\value{
  Another object of class \code{"owin"} representing the
  window after margins are added.
}
\details{
  This is a simple convenience function to add a
  margin of specified width and height on each side of a
  rectangular window. Unequal margins can also be added.
}
\seealso{
  \code{\link{trim.rectangle}},
  \code{\link{dilation}},
  \code{\link{erosion}},
  \code{\link{owin.object}}
}
\examples{
  w <- square(10)
  # add a margin of width 1 on all four sides
  square12 <- grow.rectangle(w, 1)

  # add margin of width 3 on the right side
  # and margin of height 4 on top.
  v <- grow.rectangle(w, c(0,3), c(0,4))

  # grow by 5 percent on all sides
  grow.rectangle(w, fraction=0.05)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{manip}
