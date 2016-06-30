\name{Extract.psp}
\alias{[.psp}
\title{Extract Subset of Line Segment Pattern}
\description{
  Extract a subset of a line segment pattern.
}
\usage{
  \method{[}{psp}(x, i, j, drop, ...)
}
\arguments{
  \item{x}{
    A two-dimensional line segment pattern.
    An object of class \code{"psp"}.
  }
  \item{i}{
    Subset index. Either a valid subset index in the usual \R sense,
   indicating which segments should be retained, or a window
   (an object of class \code{"owin"})
    delineating a subset of the original observation window. 
  }
  \item{j}{
   Redundant - included for backward compatibility.
  }
  \item{drop}{
    Ignored. Required for compatibility with generic function.
  }
  \item{\dots}{
    Ignored.
  }
} 
\value{
  A line segment pattern (of class \code{"psp"}).
}
\details{
  These functions extract a designated subset of a line segment pattern.

  The function \code{[.psp} is a method for \code{\link{[}} for the
  class \code{"psp"}. It extracts a designated subset of a line segment pattern,
  either by ``\emph{thinning}''
  (retaining/deleting some line segments of a line segment pattern)
  or ``\emph{trimming}'' (reducing the window of observation
  to a smaller subregion and clipping the line segments to
  this boundary) or both.
 
  The pattern will be ``thinned''
  if  \code{subset} is specified. The line segments designated by \code{subset}
  will be retained. Here \code{subset} can be a numeric vector
  of positive indices (identifying the line segments to be retained),
  a numeric vector of negative indices (identifying the line segments
  to be deleted) or a logical vector of length equal to the number
  of line segments in the line segment pattern \code{x}. In the latter case,
  the line segments for which 
  \code{subset[i]=TRUE} will be retained, and the others
  will be deleted.
 
  The pattern will be ``trimmed''
  if \code{window} is specified. This should
  be an object of class \code{\link{owin}} specifying a window of observation
  to which the line segment pattern \code{x} will be
  trimmed. Line segments of \code{x} lying inside the new
  \code{window} will be retained unchanged. Line segments lying
  partially inside the new \code{window} and partially outside it
  will be clipped so that they lie entirely inside.
 
  Both ``thinning'' and ``trimming'' can be performed together.
}
\seealso{
  \code{\link{psp.object}},
  \code{\link{owin.object}}
}
\examples{
    a <- psp(runif(20),runif(20),runif(20),runif(20), window=owin())
    plot(a)
  # thinning
    id <- sample(c(TRUE, FALSE), 20, replace=TRUE)
    b <- a[id]
    plot(b, add=TRUE, lwd=3)
 # trimming
    plot(a)
    w <- owin(c(0.1,0.7), c(0.2, 0.8))
    b <- a[w]
    plot(b, add=TRUE, col="red")
}
\author{
  \adrian
  and
  \rolf
}
\keyword{spatial}
\keyword{manip}
