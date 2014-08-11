\name{box3}
\Rdversion{1.1}
\alias{box3}
\title{
  Three-Dimensional Box
}
\description{
  Creates an object representing a three-dimensional box.
}
\usage{
box3(xrange = c(0, 1), yrange = xrange, zrange = yrange, unitname = NULL)
}
\arguments{
  \item{xrange, yrange, zrange}{
    Dimensions of the box in the \eqn{x,y,z} directions.
    Each of these arguments should be a numeric vector of length 2.
  }
  \item{unitname}{
    Optional. Name of the unit of length. See Details.
}
}
\details{
  This function creates an object representing
  a three-dimensional rectangular parallelepiped (box)
  with sides parallel to the coordinate axes.

  The object can be used to specify the domain of a three-dimensional
  point pattern (see \code{\link{pp3}}) and in various
  geometrical calculations (see \code{\link{volume.box3}},
  \code{\link{diameter.box3}}, \code{\link{eroded.volumes}}). 
  
  The optional argument \code{unitname} specifies the name
  of the unit of length. See \code{\link{unitname}}
  for valid formats.

  The function \code{\link{as.box3}} can be used to convert other kinds
  of data to this format.
}
\value{
  An object of class \code{"box3"}. There is a print method for this class.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{as.box3}}, 
  \code{\link{pp3}},
  \code{\link{volume.box3}},
  \code{\link{diameter.box3}},
  \code{\link{eroded.volumes}}.
}
\examples{
    box3()
    box3(c(0,10),c(0,10),c(0,5), unitname=c("metre","metres"))
    box3(c(-1,1))
}
\keyword{spatial}
\keyword{datagen}
