\name{rescale.owin}
\alias{rescale.owin}
\title{Convert Window to Another Unit of Length}
\description{
  Converts a window to another unit of length.
}
\usage{
  \method{rescale}{owin}(X, s, unitname)
}
\arguments{
  \item{X}{Window (object of class \code{"owin"}).}
  \item{s}{Conversion factor: the new units are \code{s} times the old units.}
  \item{unitname}{
    Optional. New name for the unit of length.
    See \code{\link{unitname}}.
  }
}
\value{
  Another window object (of class \code{"owin"}) representing the same
  window, but expressed in the new units.
}
\details{
  This is a method for the generic function \code{\link{rescale}}.

  The spatial coordinates in the window \code{X}
  (and its window) will be re-expressed
  in terms of a new unit of length that is \code{s} times the current
  unit of length given in \code{X}.
  (Thus, the coordinate values are \emph{divided} by \code{s},
  while the unit value is multiplied by \code{s}).

  The result is a window representing the \emph{same} region of space,
  but re-expressed in a different unit.

  If \code{s} is missing, then the coordinates will be re-expressed
  in \sQuote{native} units; for example if the current unit is
  equal to 0.1 metres, then the coordinates will be re-expressed in metres.
}
\section{Note}{
  The result of this operation is equivalent to the original window.
  If you want to actually change the coordinates by
  a linear transformation, producing a window that is larger or smaller
  than the original one, use \code{\link{affine}}.
}
\seealso{
  \code{\link{unitname}},
  \code{\link{rescale}},
  \code{\link{rescale.owin}},
  \code{\link{affine}},
  \code{\link{rotate}},
  \code{\link{shift}}
}
\examples{
  data(swedishpines)
  W <- swedishpines$window
  W
# coordinates are in decimetres (0.1 metre)
# convert to metres:
  rescale(W, 10)
# or equivalently
  rescale(W)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
