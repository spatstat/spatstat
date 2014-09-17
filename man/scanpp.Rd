\name{scanpp}
\alias{scanpp}
\title{Read Point Pattern From Data File} 
\description{
Reads a point pattern dataset from a text file.
}
\usage{
   scanpp(filename, window, header=TRUE, dir="", factor.marks=NULL, ...) 
}
\arguments{
  \item{filename}{
    String name of the file containing
    the coordinates of the points in the point pattern,
    and their marks if any.
    }
    \item{window}{
      Window for the point pattern. An object of class \code{"owin"}.
    }
    \item{header}{
      Logical flag indicating whether the first line of the
      file contains headings for the columns.
      Passed to \code{\link[utils]{read.table}}.
    }
    \item{dir}{
      String containing the path name of the directory
      in which \code{filename} is to be found.
      Default is the current directory.
    }
    \item{factor.marks}{
      Logical vector (or NULL) indicating whether marks are to be
      interpreted as factors. Defaults to \code{NULL} which means that
      strings will be interpreted as factors while numeric variables
      will not. See details.
    }
    \item{\dots}{
      Ignored.
    }
}
\value{
  A point pattern (an object of class \code{"ppp"},
  see \code{\link{ppp.object}}).
}
\details{
  This simple function reads a point pattern dataset from a file
  containing the cartesian coordinates of its points,
  and optionally the mark values for these points.

  The file identified by \code{filename} in directory \code{dir}
  should be a text file that can be read using \code{\link[utils]{read.table}}.
  Thus, each line of the file (except possibly the first line)
  contains data for one point in the
  point pattern. Data are arranged in columns. There should be either
  two columns (for an unmarked point pattern) or more columns (for a
  marked point pattern).

  If \code{header=FALSE} then the first two columns of data
  will be interpreted as the \eqn{x} and \eqn{y} coordinates
  of points. Remaining columns, if present, will be interpreted as
  containing the marks for these points.

  If \code{header=TRUE} then the first line of the file should contain
  string names for each of the columns of data. If there are columns
  named \code{x} and \code{y} then these will be taken as the
  cartesian coordinates, and any remaining columns will be taken as
  the marks. If there are no columns named \code{x} and \code{y}
  then the first and second columns will be taken as the cartesian
  coordinates.

  If a logical vector is provided for \code{factor.marks} the length
  should equal the number of mark columns (a shorter \code{factor.marks}
  is recycled to this length). This vector is then used to determine
  which mark columns should be interpreted as factors. Note: Strings will
  not be interpreted as factors if the corresponding entry in
  \code{factor.marks} is \code{FALSE}.

  Note that there is intentionally no default for \code{window}.
  The window of observation should be specified.
  If you really need to estimate the window, use the
  Ripley-Rasson estimator \code{\link{ripras}}.
}
\seealso{
  \code{\link{ppp.object}},
  \code{\link{ppp}},
  \code{\link{as.ppp}},
  \code{\link{ripras}}
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{IO}
 
 
