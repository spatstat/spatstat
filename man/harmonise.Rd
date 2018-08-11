\name{harmonise}
\alias{harmonise}
\alias{harmonize}
\title{Make Objects Compatible}
\description{
  Converts several objects of the same class to a common format
  so that they can be combined or compared.
}
\usage{
harmonise(\dots)
harmonize(\dots)
}
\arguments{
  \item{\dots}{
    Any number of objects of the same class.
  } 
}
\details{
  This generic command takes any number of objects of the same
  class, and \emph{attempts} to make them compatible in the sense
  of \code{\link{compatible}} so that they can be combined or compared.

  There are methods for the classes \code{"fv"}
  (\code{\link{harmonise.fv}})
  and \code{"im"} (\code{\link{harmonise.im}}).
  
  All arguments \code{\dots} must be objects of the same class.
  The result will be a list, of length equal to the number of
  arguments \code{\dots}, containing new versions of each of these
  objects, converted to a common format.
  If the arguments were named (\code{name=value}) then the return value
  also carries these names.
}
\value{
  A list, of length equal to the number of arguments \code{\dots},
  whose entries are objects of the same class.
  If the arguments were named (\code{name=value}) then the return value
  also carries these names.
}
\author{\adrian
  
  ,
  \rolf
  
  and \ege.
}
\seealso{
  \code{\link{compatible}},
  \code{\link{harmonise.fv}},
  \code{\link{harmonise.im}}
}
\keyword{spatial}
\keyword{manip}
