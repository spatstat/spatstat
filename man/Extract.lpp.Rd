\name{Extract.lpp}
\alias{[.lpp}
\title{Extract Subset of Point Pattern on Linear Network}
\description{
  Extract a subset of a point pattern on a linear network.
}
\usage{
  \method{[}{lpp}(x, i, j, drop=FALSE, \dots, snip=TRUE)
}
\arguments{
  \item{x}{
    A point pattern on a linear network (object of class \code{"lpp"}).
  }
  \item{i}{
   Subset index. A valid subset index in the usual \R sense,
   indicating which points should be retained.
  }
  \item{j}{
    Spatial window
   (object of class \code{"owin"}) delineating the region that should
   be retained.
  }
  \item{drop}{
    Logical value indicating whether to remove unused levels
    of the marks, if the marks are a factor.
  }
  \item{snip}{
    Logical. If \code{TRUE} (the default), segments of the network
    which cross the boundary of the window \code{j} will be cut by the boundary.
    If \code{FALSE}, these segments will be deleted.
  }
  \item{\dots}{
    Ignored. 
  }
}
\value{
  A point pattern on a linear network (of class \code{"lpp"}).
}
\details{
  This function extracts a designated subset of a point pattern
  on a linear network.

  The function \code{[.lpp} is a method for \code{\link{[}} for the
  class \code{"lpp"}. It extracts a designated subset of a point pattern.
  The argument \code{i} should be a subset index in the usual \R sense:
  either a numeric vector
  of positive indices (identifying the points to be retained),
  a numeric vector of negative indices (identifying the points
  to be deleted) or a logical vector of length equal to the number of
  points in the point pattern \code{x}. In the latter case, 
  the points \code{(x$x[i], x$y[i])} for which 
  \code{subset[i]=TRUE} will be retained, and the others
  will be deleted.

  The argument \code{j}, if present, should be a spatial window.
  The pattern inside the region will be retained.
  \emph{Line segments that cross the boundary of the window
  are deleted} in the current implementation.

  The argument \code{drop} determines whether to remove
  unused levels of a factor, if the point pattern is multitype
  (i.e. the marks are a factor) or if the marks are a data frame or hyperframe
  in which some of the columns are factors.

  The argument \code{snip} specifies what to do with segments of
  the network which cross the boundary of the window \code{j}. 
  If \code{snip=FALSE}, such segments are simply deleted.  
  If \code{snip=TRUE} (the default), such segments are cut into pieces by the
  boundary of \code{j}, and those pieces which lie inside the window
  \code{ji} are included in the resulting network. 

  Use \code{\link{unmark}} to remove all the marks in a marked point
  pattern, and \code{\link{subset.lpp}} to remove only some columns of marks. 
}
\seealso{
  \code{\link{lpp}},
  \code{\link{subset.lpp}}
}
\examples{
  # Chicago crimes data - remove cases of assault
  chicago[marks(chicago) != "assault"]
  # equivalent to subset(chicago, select=-assault)

  # spatial window subset
  B <- owin(c(350, 700), c(600, 1000))
  plot(chicago)
  plot(B, add=TRUE, lty=2, border="red", lwd=3)
  op <- par(mfrow=c(1,2), mar=0.6+c(0,0,1,0))
  plot(B, main="chicago[B, snip=FALSE]", lty=3, border="red")
  plot(chicago[, B, snip=FALSE], add=TRUE)
  plot(B, main="chicago[B, snip=TRUE]", lty=3, border="red")
  plot(chicago[, B, snip=TRUE], add=TRUE)
  par(op)
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  and
  Rolf Turner \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
