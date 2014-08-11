\name{rescale.ppp}
\alias{rescale.ppp}
\title{Convert Point Pattern to Another Unit of Length}
\description{
  Converts a point pattern dataset to 
  another unit of length.
}
\usage{
 \method{rescale}{ppp}(X, s)
}
\arguments{
  \item{X}{Point pattern (object of class \code{"ppp"}).}
  \item{s}{Conversion factor: the new units are \code{s} times the old units.}
}
\value{
  Another point pattern (of class \code{"ppp"}),
  representing the same data, but expressed in the new units.
}
\details{
  This is a method for the generic function \code{\link{rescale}}.

  The spatial coordinates in the point pattern \code{X}
  (and its window) will be re-expressed
  in terms of a new unit of length that is \code{s} times the current
  unit of length given in \code{X}.
  (Thus, the coordinate values are \emph{divided} by \code{s},
  while the unit value is multiplied by \code{s}).

  The result is a point pattern representing the \emph{same} data
  but re-expressed in a different unit.

  Mark values are unchanged.

  If \code{s} is missing, then the coordinates will be re-expressed
  in \sQuote{native} units; for example if the current unit is
  equal to 0.1 metres, then the coordinates will be re-expressed in metres.
}
\section{Note}{
  The result of this operation is equivalent to the original point pattern.
  If you want to actually change the coordinates by
  a linear transformation, producing a point pattern that is not
  equivalent to the original one, use \code{\link{affine}}.
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
# Bramble Canes data: 1 unit = 9 metres
  data(bramblecanes)
# convert to metres
  bram <- rescale(bramblecanes, 1/9)
# or equivalently
  bram <- rescale(bramblecanes)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
