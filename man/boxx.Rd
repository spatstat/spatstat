\name{boxx}
\Rdversion{1.1}
\alias{boxx}
\title{
  Multi-Dimensional Box
}
\description{
  Creates an object representing a multi-dimensional box.
}
\usage{
boxx(..., unitname = NULL)
}
\arguments{
  \item{\dots}{
    Dimensions of the box. Vectors of length 2.
  }
  \item{unitname}{
    Optional. Name of the unit of length. See Details.
  }
}
\details{
  This function creates an object representing
  a multi-dimensional rectangular parallelepiped (box)
  with sides parallel to the coordinate axes.

  The object can be used to specify the domain of a multi-dimensional
  point pattern (see \code{\link{ppx}}) and in various
  geometrical calculations (see \code{\link{volume.boxx}},
  \code{\link{diameter.boxx}}, \code{\link{eroded.volumes}}). 
  
  The optional argument \code{unitname} specifies the name
  of the unit of length. See \code{\link{unitname}}
  for valid formats.
}
\value{
  An object of class \code{"boxx"}. There is a print method for this class.
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\seealso{
  \code{\link{ppx}},
  \code{\link{volume.boxx}},
  \code{\link{diameter.boxx}},
  \code{\link{eroded.volumes.boxx}}.
}
\examples{
    boxx(c(0,10),c(0,10),c(0,5),c(0,1), unitname=c("metre","metres"))
}
\keyword{spatial}
\keyword{datagen}
