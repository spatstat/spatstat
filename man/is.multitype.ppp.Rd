\name{is.multitype.ppp}
\alias{is.multitype.ppp}
\alias{is.multitype.lpp}
\title{Test Whether A Point Pattern is Multitype}
\description{
  Tests whether a point pattern has ``marks''
  attached to the points which classify the points into several types.
}
\usage{
  \method{is.multitype}{ppp}(X, na.action="warn", \dots) 

  \method{is.multitype}{lpp}(X, na.action="warn", \dots) 
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"} or \code{"lpp"})
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
  Logical value, equal to \code{TRUE} if \code{X} is a multitype point pattern.
}
\details{
  ``Marks'' are observations attached to each point of a point pattern.
  For example the \code{\link[spatstat.data]{longleaf}} dataset
  contains the locations of trees, each tree being marked by its diameter;
  the \code{\link[spatstat.data]{amacrine}} dataset gives the locations of cells
  of two types (on/off) and the type of cell may be regarded as a mark attached
  to the location of the cell.

  This function tests whether the point pattern \code{X}
  contains or involves marked points, \bold{and} that the
  marks are a factor.
  It is a method for the generic function \code{\link{is.multitype}}.
  
  For example, the \code{\link[spatstat.data]{amacrine}}
  dataset is multitype (there are two types of cells, on and off),
  but the \code{\link[spatstat.data]{longleaf}} dataset is \emph{not} multitype
  (the marks are real numbers).

  The argument \code{na.action} determines what action will be taken
  if the point pattern has a vector of marks but some or all of the
  marks are \code{NA}. Options are   \code{"fatal"} to cause a fatal
  error; \code{"warn"} to issue a warning and then return \code{TRUE};
  and \code{"ignore"} to take no action except returning \code{TRUE}.
}
\seealso{
  \code{\link{is.multitype}},
  \code{\link{is.multitype.ppm}}
}
\examples{
   is.multitype(cells)  #FALSE - no marks
   is.multitype(longleaf) #FALSE - real valued marks
   is.multitype(amacrine) #TRUE
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{manip}
