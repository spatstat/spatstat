\name{is.marked.ppp}
\alias{is.marked.ppp}
\title{Test Whether A Point Pattern is Marked}
\description{
  Tests whether a point pattern has ``marks''
  attached to the points.
}
\usage{
  \method{is.marked}{ppp}(X, na.action="warn", \dots) 
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"})
  }
  \item{na.action}{
    String indicating what to do if \code{NA} values are
    encountered amongst the marks.
    Options are \code{"warn"}, \code{"fatal"} and \code{"ignore"}.
  }
  \item{\dots}{
    Ignored.
  }
}
\value{
  Logical value, equal to \code{TRUE} if \code{X} is a marked point pattern.
}
\details{
  ``Marks'' are observations attached to each point of a point pattern.
  For example the \code{\link[spatstat.data]{longleaf}} dataset
  contains the locations of trees, each tree being marked by its diameter;
  the \code{\link[spatstat.data]{amacrine}} dataset gives the locations of cells
  of two types (on/off) and the type of cell may be regarded as a mark attached
  to the location of the cell.

  This function tests whether the point pattern \code{X}
  contains or involves marked points.
  It is a method for the generic function \code{\link{is.marked}}.

  The argument \code{na.action} determines what action will be taken
  if the point pattern has a vector of marks but some or all of the
  marks are \code{NA}. Options are   \code{"fatal"} to cause a fatal
  error; \code{"warn"} to issue a warning and then return \code{TRUE};
  and \code{"ignore"} to take no action except returning \code{TRUE}.
}
\seealso{
  \code{\link{is.marked}},
  \code{\link{is.marked.ppm}}
}
\examples{
   data(cells)
   is.marked(cells)  #FALSE
   data(longleaf)
   is.marked(longleaf) #TRUE
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{manip}
