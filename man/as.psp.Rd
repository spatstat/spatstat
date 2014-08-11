\name{as.psp}
\alias{as.psp}
\alias{as.psp.psp}
\alias{as.psp.data.frame}
\alias{as.psp.matrix}
\alias{as.psp.owin}
\alias{as.psp.default}
\title{Convert Data To Class psp}
\description{
  Tries to coerce any reasonable kind of data object to a line segment pattern
  (an object of class \code{"psp"})
  for use by the \pkg{spatstat} package.
}
\usage{
  as.psp(x, \dots, from=NULL, to=NULL)

  \method{as.psp}{psp}(x, \dots, check=FALSE, fatal=TRUE)

  \method{as.psp}{data.frame}(x, \dots, window=NULL, marks=NULL,
      check=spatstat.options("checksegments"), fatal=TRUE)

  \method{as.psp}{matrix}(x, \dots, window=NULL, marks=NULL,
       check=spatstat.options("checksegments"), fatal=TRUE)

  \method{as.psp}{owin}(x, \dots,
       check=spatstat.options("checksegments"), fatal=TRUE)

  \method{as.psp}{default}(x, \dots, window=NULL, marks=NULL,
       check=spatstat.options("checksegments"), fatal=TRUE)
}
\arguments{
  \item{x}{Data which will be converted into a line segment pattern}
  \item{window}{Data which define a window for the pattern
    when \code{x} does not contain a window}
  \item{\dots}{Ignored.}
  \item{marks}{(Optional) vector or data frame of marks for the pattern}
  \item{check}{
    Logical value indicating whether to check the validity of the data,
    e.g. to check that the line segments lie inside the window.
  }
  \item{fatal}{Logical value. See Details.}
  \item{from,to}{Point patterns (object of class \code{"ppp"})
    containing the first and second endpoints (respectively) of each
    segment. Incompatible with \code{x}.
  }
}
\value{
  An object of class \code{"psp"} (see \code{\link{psp.object}})
  describing the line segment pattern and its window of observation.
  The value \code{NULL} may also be returned; see Details.
}
\details{
  Converts the dataset \code{x} to a line segment pattern
  (an object of class \code{"psp"}; see \code{\link{psp.object}} for
  an overview).

  This function is normally used to convert an existing line segment pattern
  dataset, stored in another format, to the \code{"psp"} format.
  To create a new point pattern from raw data such as \eqn{x,y}
  coordinates, it is normally easier to use the creator function
  \code{\link{psp}}.

  The dataset \code{x} may be:
  \itemize{
    \item
    an object of class \code{"psp"}  
    \item
    a data frame with at least 4 columns
    \item a structure (list) with elements named  \code{x0, y0,
    x1, y1} or elements named \code{xmid, ymid, length, angle}
    and possibly a fifth element named \code{marks}
    \item
    an object of class \code{"owin"} representing a spatial window;
    it must be of type \code{"rectangle"} or \code{"polygonal"}.
    The boundary edges of the window will be extracted as a line
    segment pattern.
  }
  If \code{x} is a data frame the interpretation of its columns is
  as follows:
  \itemize{
    \item
    If there are columns named \code{x0, y0, x1, y1} then these
    will be interpreted as the coordinates of the endpoints of
    the segments and used to form the \code{ends} component of
    the \code{psp} object to be returned.
    \item
    If there are columns named \code{xmid, ymid, length, angle}
    then these will be interpreted as the coordinates of the segment
    midpoints, the lengths of the segments, and the orientations
    of the segments in radians and used to form the \code{ends}
    component of the \code{psp} object to be returned.
    \item
    If there is a column named \code{marks} then this will
    be interpreted as the marks of the pattern provided that
    the argument \code{marks} of this function is \code{NULL}.
    If argument \code{marks} is not \code{NULL} then the value
    of this argument is taken to be the marks of the pattern and
    the column named \code{marks} is ignored (with a warning).
    In either case the column named marks is deleted and omitted
    from further consideration.
    \item
    If there is no column named \code{marks} and if the \code{marks}
    argument of this function is \code{NULL}, and if after
    interpreting 4 columns of \code{x} as determining the \code{ends}
    component of the \code{psp} object to be returned, there remain
    other columns of \code{x}, then these remaining columns will be
    taken to form a data frame of marks for the \code{psp} object
    to be returned.
  }

  If \code{x} is a structure (list) with elements named \code{x0,
  y0, x1, y1, marks} or \code{xmid, ymid, length, angle, marks},
  then the element named \code{marks} will be interpreted as the
  marks of the pattern provide that the argument \code{marks} of
  this function is \code{NULL}.  If this argument is non-\code{NULL}
  then it is interpreted as the marks of the pattern and the element
  \code{marks} of \code{x} is ignored --- with a warning.

  Alternatively, you may specify two point patterns
  \code{from} and \code{to} containing the first and second endpoints
  of the line segments.
  
  The argument \code{window} is converted to a window object
  by the function \code{\link{as.owin}}.

  The argument \code{fatal} indicates what to do when
  the data cannot be converted to a line segment pattern.
  If \code{fatal=TRUE}, a fatal error
  will be generated; if \code{fatal=FALSE}, the
  value \code{NULL} is returned.

  The function \code{as.psp} is generic, with methods for the
  classes \code{"psp"}, \code{"data.frame"}, \code{"matrix"}
  and a default method.
  
  Point pattern datasets can also be created by the function
  \code{\link{psp}}.
}
\section{Warnings}{
If only a proper subset of the names \code{x0,y0,x1,y1} or
\code{xmid,ymid,length,angle} appear amongst the names of the
columns of \code{x} where \code{x} is a data frame, then these
special names are ignored.

For example if the names of the columns were
\code{xmid,ymid,length,degrees}, then these columns would be
interpreted as if the represented \code{x0,y0,x1,y1} in that order.

Whether it gets used or not, column named \code{marks} is
\emph{always} removed from \code{x} before any attempt to form the
\code{ends} component of the \code{psp} object that is returned.
}
\seealso{
  \code{\link{psp}},
  \code{\link{psp.object}},
  \code{\link{as.owin}},
  \code{\link{owin.object}}
}
\examples{
   mat <- matrix(runif(40), ncol=4)
   mx <- data.frame(v1=sample(1:4,10,TRUE),
                    v2=factor(sample(letters[1:4],10,TRUE),levels=letters[1:4]))
   a <- as.psp(mat, window=owin(),marks=mx)
   mat <- cbind(as.data.frame(mat),mx)
   b <- as.psp(mat, window=owin()) # a and b are identical.
   stuff <- list(xmid=runif(10),
                 ymid=runif(10),
                 length=rep(0.1, 10),
                 angle=runif(10, 0, 2 * pi))
   a <- as.psp(stuff, window=owin())
   b <- as.psp(from=runifpoint(10), to=runifpoint(10))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
