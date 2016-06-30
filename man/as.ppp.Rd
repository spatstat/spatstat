\name{as.ppp}
\alias{as.ppp}
\alias{as.ppp.ppp}
\alias{as.ppp.psp}
\alias{as.ppp.quad}
\alias{as.ppp.matrix}
\alias{as.ppp.data.frame}
\alias{as.ppp.influence.ppm}
\alias{as.ppp.default}
\title{Convert Data To Class ppp}
\description{
  Tries to coerce any reasonable kind of data to a spatial point pattern
  (an object of class \code{"ppp"})
  for use by the \pkg{spatstat} package).
}
\usage{
  as.ppp(X, \dots, fatal=TRUE)

  \method{as.ppp}{ppp}(X, \dots, fatal=TRUE)

  \method{as.ppp}{psp}(X, \dots, fatal=TRUE)

  \method{as.ppp}{quad}(X, \dots, fatal=TRUE)

  \method{as.ppp}{matrix}(X, W=NULL, \dots, fatal=TRUE)

  \method{as.ppp}{data.frame}(X, W=NULL, \dots, fatal=TRUE)

  \method{as.ppp}{influence.ppm}(X, \dots)

  \method{as.ppp}{default}(X, W=NULL, \dots, fatal=TRUE)
}
\arguments{
  \item{X}{Data which will be converted into a point pattern}
  \item{W}{
    Data which define a window for the pattern,
    when \code{X} does not contain a window.
    (Ignored if \code{X} contains window information.)
  }
  \item{\dots}{Ignored.}
  \item{fatal}{
    Logical value specifying what to do if the
    data cannot be converted.
    See Details.
  }
}
\value{
  An object of class \code{"ppp"} (see \code{\link{ppp.object}})
  describing the point pattern and its window of observation.
  The value \code{NULL} may also be returned; see Details.
}
\details{
  Converts the dataset \code{X} to a point pattern
  (an object of class \code{"ppp"}; see \code{\link{ppp.object}} for
  an overview).

  This function is normally used to convert an existing point pattern
  dataset, stored in another format, to the \code{"ppp"} format.
  To create a new point pattern from raw data such as \eqn{x,y}
  coordinates, it is normally easier to use the creator function
  \code{\link{ppp}}.

  The function \code{as.ppp} is generic, with methods for the
  classes \code{"ppp"}, \code{"psp"}, \code{"quad"}, \code{"matrix"},
  \code{"data.frame"} and a default method.

  The dataset \code{X} may be:
  \itemize{
    \item
    an object of class \code{"ppp"}  
    \item
    an object of class \code{"psp"}
    \item
    a point pattern object created by the \pkg{spatial} library
    \item
    an object of class \code{"quad"} representing a quadrature scheme
    (see \code{\link{quad.object}})
    \item
    a matrix or data frame with at least two columns
    \item
    a structure with entries \code{x}, \code{y} which are numeric vectors
    of equal length
    \item
    a numeric vector of length 2, interpreted as the coordinates
    of a single point.
  }
  In the last three cases, we need the second argument \code{W}
  which is converted to a window object
  by the function \code{\link{as.owin}}.
  In the first four cases, \code{W} will be ignored.

  If \code{X} is a line segment pattern (an object of class \code{psp})
  the point pattern returned consists of the endpoints of the segments.
  If \code{X} is marked then the point pattern returned will also be
  marked, the mark associated with a point being the mark of the segment
  of which that point was an endpoint.

  If \code{X} is a matrix or data frame, the first and second columns will
  be interpreted as the \eqn{x} and \eqn{y} coordinates respectively.
  Any additional columns will be interpreted as marks.
  
  The argument \code{fatal} indicates what to do when 
  \code{W} is missing and \code{X} contains no
  information about the window. If \code{fatal=TRUE}, a fatal error
  will be generated; if \code{fatal=FALSE}, the
  value \code{NULL} is returned.

  In the \pkg{spatial} library, a point pattern is represented
  in either of the following formats:
  \itemize{
    \item
    (in \pkg{spatial} versions 1 to 6)
    a structure with entries \code{x}, \code{y}
    \code{xl}, \code{xu}, \code{yl}, \code{yu}
    \item
    (in \pkg{spatial} version 7)
    a structure with entries
    \code{x}, \code{y} and \code{area},
    where \code{area} is a structure with entries
    \code{xl}, \code{xu}, \code{yl}, \code{yu}
  }
  where \code{x} and \code{y} are vectors of equal length
  giving the point coordinates, and \code{xl}, \code{xu}, \code{yl},
  \code{yu} are numbers giving the dimensions of a rectangular window.

  Point pattern datasets can also be created by the function
  \code{\link{ppp}}.
}
\seealso{
  \code{\link{ppp}}, \code{\link{ppp.object}}, \code{\link{as.owin}},
  \code{\link{owin.object}}
}
\examples{
 xy <- matrix(runif(40), ncol=2)
 pp <- as.ppp(xy, c(0,1,0,1))

 # Venables-Ripley format
 # check for 'spatial' package
 spatialpath <- system.file(package="spatial")
 if(nchar(spatialpath) > 0) {
   require(spatial)
   towns <- ppinit("towns.dat")
   pp <- as.ppp(towns) # converted to our format
   detach(package:spatial)
 }

 xyzt <- matrix(runif(40), ncol=4)
 Z <- as.ppp(xyzt, square(1))
}
\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{manip}
