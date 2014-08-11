\name{as.lpp}
\Rdversion{1.1}
\alias{as.lpp}
\title{
  Convert Data to a Point Pattern on a Linear Network
}
\description{
  Convert various kinds of data to a point pattern on a linear network.
}
\usage{
  as.lpp(x, y=NULL, seg=NULL, tp=NULL, \dots,
         marks=NULL, L=NULL, check=FALSE)
}
\arguments{
  \item{x,y}{
    Vectors of cartesian coordinates, or any data
    acceptable to \code{\link[grDevices]{xy.coords}}.
    Alternatively \code{x} can be a point pattern
    on a linear network (object of class \code{"lpp"})
    or a planar point pattern (object of class \code{"ppp"}).
  }
  \item{seg,tp}{
    Optional local coordinates. Vectors of the same length
    as \code{x,y}. See Details.
  }
  \item{\dots}{Ignored.}
  \item{marks}{
    Optional marks for the point pattern.
    A vector or factor with one entry for each point,
    or a data frame or hyperframe with one row for each point.
  }
  \item{L}{
    Linear network (object of class \code{"linnet"})
    on which the points lie.
  }
  \item{check}{
    Logical. Whether to check the validity of the spatial coordinates.
  }
}
\details{
  This function converts data in various formats into a point pattern
  on a linear network  (object of class \code{"lpp"}).
  
  The possible formats are:
  \itemize{
    \item
    \code{x} is already a point pattern on a linear network
    (object of class \code{"lpp"}).
    Then \code{x} is returned unchanged.
    \item
    \code{x} is a planar point pattern (object of class \code{"ppp"}).
    Then \code{x} is converted to a point pattern on the linear network
    \code{L} using \code{\link{lpp}}.
    \item
    \code{x,y,seg,tp} are vectors of equal length.
    These specify that the \code{i}th point has Cartesian coordinates
    \code{(x[i],y[i])}, and lies on segment number \code{seg[i]} of the
    network \code{L}, at a fractional position \code{tp[i]} along that
    segment (with \code{tp=0} representing one endpoint and
    \code{tp=1} the other endpoint of the segment).
    \item
    \code{seg,tp} are \code{NULL}, and \code{x,y} are data in a format
    acceptable to \code{\link[grDevices]{xy.coords}} specifying the
    Cartesian coordinates.
  }
}
\value{
  A point pattern
  on a linear network (object of class \code{"lpp"}).
}
\seealso{
  \code{\link{lpp}}.
}
\examples{
   A <- as.psp(simplenet)
   X <- runifpointOnLines(10, A)
   is.ppp(X)
   Y <- as.lpp(X, L=simplenet)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
