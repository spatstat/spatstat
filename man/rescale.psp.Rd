\name{rescale.psp}
\alias{rescale.psp}
\title{Convert Line Segment Pattern to Another Unit of Length}
\description{
  Converts a line segment pattern dataset to 
  another unit of length.
}
\usage{
 \method{rescale}{psp}(X, s, unitname)
}
\arguments{
  \item{X}{Line segment pattern (object of class \code{"psp"}).}
  \item{s}{Conversion factor: the new units are \code{s} times the old units.}
  \item{unitname}{
    Optional. New name for the unit of length.
    See \code{\link{unitname}}.
  }
}
\value{
  Another line segment pattern (of class \code{"psp"}),
  representing the same data, but expressed in the new units.
}
\details{
  This is a method for the generic function \code{\link{rescale}}.

  The spatial coordinates in the line segment pattern \code{X}
  (and its window) will be re-expressed
  in terms of a new unit of length that is \code{s} times the current
  unit of length given in \code{X}.
  (Thus, the coordinate values are \emph{divided} by \code{s},
  while the unit value is multiplied by \code{s}).

  The result is a line segment pattern representing the \emph{same} data
  but re-expressed in a different unit.

  Mark values are unchanged.

  If \code{s} is missing, then the coordinates will be re-expressed
  in \sQuote{native} units; for example if the current unit is
  equal to 0.1 metres, then the coordinates will be re-expressed in metres.
}
\section{Note}{
  The result of this operation is equivalent to the original segment pattern.
  If you want to actually change the coordinates by
  a linear transformation, producing a segment pattern that is not
  equivalent to the original one, use \code{\link{affine}}.
}
\seealso{
  \code{\link{units}},
  \code{\link{affine}},
  \code{\link{rotate}},
  \code{\link{shift}}
}
\examples{
   data(copper)
   X <- copper$Lines
   X
   # data are in km
   # convert to metres
   rescale(X, 1/1000)

   # convert data and rename unit
   rescale(X, 1/1000, c("metre", "metres"))
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{math}
