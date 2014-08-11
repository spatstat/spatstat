\name{scanpp}
\alias{scanpp}
\title{Read Point Pattern From Data File} 
\description{
Reads a point pattern dataset from a text file.
}
\usage{
   scanpp(filename, window, header=TRUE, dir="", multitype=FALSE) 
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
      Passed to \code{\link{read.table}}.
    }
    \item{dir}{
      String containing the path name of the directory
      in which \code{filename} is to be found.
      Default is the current directory.
    }
    \item{multitype}{
      Logical flag indicating whether marks are to be
      interpreted as a factor (\code{multitype = TRUE})
      or as numerical values (\code{multitype = FALSE}).
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
  should be a text file that can be read using \code{\link{read.table}}.
  Thus, each line of the file (except possibly the first line)
  contains data for one point in the
  point pattern. Data are arranged in columns. There should be either
  two columns (for an unmarked point pattern) or three columns (for a
  marked point pattern).

  If \code{header=FALSE} then the first two columns of data
  will be interpreted as the \eqn{x} and \eqn{y} coordinates
  of points. A third column, if present, will be interpreted as
  containing the marks for these points. The marks will be
  converted to a factor if \code{multitype = TRUE}.

  If \code{header=TRUE} then the first line of the file should contain
  string names for each of the columns of data. If there are columns
  named \code{x} and \code{y} then these will be taken as the
  cartesian coordinates, and any remaining column will be taken as
  the marks. If there are no columns named \code{x} and \code{y}
  then the first and second columns will be taken as the cartesian
  coordinates.

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
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{IO}
 
 
