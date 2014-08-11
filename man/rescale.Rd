\name{rescale}
\alias{rescale}
\title{Convert dataset to another unit of length}
\description{
  Converts between different units of length
  in a spatial dataset, such as a point pattern or a window. 
}
\usage{
  rescale(X, s, unitname)
}
\arguments{
  \item{X}{Any suitable dataset representing a two-dimensional
    object, such as a point pattern (object of class \code{"ppp"}),
    or a window (object of class \code{"owin"}).}
  \item{s}{Conversion factor: the new units are \code{s} times the old
    units.}
  \item{unitname}{
    Optional. New name for the unit of length.
    See \code{\link{unitname}}.
  }
}
\value{
  Another object of the same type, representing the same
  data, but expressed in the new units.
}
\details{
  This is generic. Methods are provided for many spatial objects.

  The spatial coordinates in the dataset \code{X} will be re-expressed
  in terms of a new unit of length that is \code{s} times the current
  unit of length given in \code{X}.

  If \code{unitname} is given, it will be taken as the new name of the unit
  of length. It should be a valid name for the
  unit of length, as described in the help for \code{\link{unitname}}.

  For example if \code{X} is a dataset giving coordinates in metres,
  then \code{rescale(X,1000)}
  will divide the coordinate values by 1000 to obtain coordinates in
  kilometres, and the unit name will be changed from \code{"metres"}
  to \code{"1000 metres"}.
  Similarly \code{rescale(X, 1000, "km")} will divide the coordinate
  values by 1000 to obtain coordinates in
  kilometres, and the unit name will be changed to \code{"km"}.
}
\section{Note}{
  The result of this operation is equivalent to the original dataset.
  If you want to actually change the coordinates by
  a linear transformation, producing a dataset that is not equivalent
  to the original one, use \code{\link{affine}}.
}
\seealso{
Available methods:
\code{\link{rescale.im}},
\code{\link{rescale.layered}},
\code{\link{rescale.linnet}},
\code{\link{rescale.lpp}},
\code{\link{rescale.owin}},
\code{\link{rescale.ppp}},
\code{\link{rescale.psp}}
and \code{\link{rescale.units}}.

Other generics: 
\code{\link{unitname}},
  \code{\link{affine}},
  \code{\link{rotate}},
  \code{\link{shift}}.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
