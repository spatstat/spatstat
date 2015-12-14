\name{Extract.ppp}
\alias{[.ppp}
\alias{[<-.ppp}
\title{Extract or Replace Subset of Point Pattern}
\description{
  Extract or replace a subset of a point pattern.
  Extraction of a subset has the effect of thinning the 
  points and/or trimming the window.
}
\usage{
  \method{[}{ppp}(x, i, j, drop=FALSE, \dots, clip=FALSE)
  \method{[}{ppp}(x, i, j) <- value
}
\arguments{
  \item{x}{
    A two-dimensional point pattern.
    An object of class \code{"ppp"}.
  }
  \item{i}{
   Subset index. Either a valid subset index in the usual \R sense,
   indicating which points should be retained, or a window
   (an object of class \code{"owin"}) 
   delineating a subset of the original observation window,
   or a pixel image with logical values defining a subset of the
   original observation window.
  }
  \item{value}{
    Replacement value for the subset. A point pattern.
  }
  \item{j}{
    Redundant. Included for backward compatibility.
  }
  \item{drop}{
    Logical value indicating whether to remove unused levels
    of the marks, if the marks are a factor.
  }
  \item{clip}{
    Logical value indicating how to form the window of the resulting
    point pattern, when \code{i} is a window. 
    If \code{clip=FALSE} (the default), the result has window
    equal to \code{i}. If \code{clip=TRUE}, the resulting window
    is the intersection between the window of \code{x} and the
    window \code{i}.
  }
  \item{\dots}{
    Ignored. This argument is required for compatibility
    with the generic function.
  }
}
\value{
  A point pattern (of class \code{"ppp"}).
}
\details{
  These functions extract a designated subset of a point pattern,
  or replace the designated subset with another point pattern.

  The function \code{[.ppp} is a method for \code{\link{[}} for the
  class \code{"ppp"}. It extracts a designated subset of a point pattern,
  either by ``\emph{thinning}''
  (retaining/deleting some points of a point pattern)
  or ``\emph{trimming}'' (reducing the window of observation
  to a smaller subregion and retaining only
  those points which lie in the subregion) or both.

  The pattern will be ``thinned''
  if \code{i} is a subset index in the usual \R sense:
  either a numeric vector
  of positive indices (identifying the points to be retained),
  a numeric vector of negative indices (identifying the points
  to be deleted) or a logical vector of length equal to the number of
  points in the point pattern \code{x}. In the latter case, 
  the points \code{(x$x[i], x$y[i])} for which 
  \code{subset[i]=TRUE} will be retained, and the others
  will be deleted.
 
  The pattern will be ``trimmed''
  if \code{i} is an object of class 
  \code{"owin"} specifying a window of observation.
  The points of \code{x} lying inside the new
  window \code{i} will be retained. Alternatively \code{i} may be a
  pixel image (object of class \code{"im"}) with logical values;
  the pixels with the value \code{TRUE} will be interpreted as a window.

  The argument \code{drop} determines whether to remove
  unused levels of a factor, if the point pattern is multitype
  (i.e. the marks are a factor) or if the marks are a data frame
  in which some of the columns are factors.

  The function \code{[<-.ppp} is a method for \code{\link{[<-}} for the
  class \code{"ppp"}. It replaces the designated
  subset with the point pattern \code{value}.
  The subset of \code{x} to be replaced is designated by
  the argument \code{i} as above.

  The replacement point pattern \code{value} must lie inside the
  window of the original pattern \code{x}.
  The ordering of points in \code{x} will be preserved
  if the replacement pattern \code{value} has the same number of points
  as the subset to be replaced.  Otherwise the ordering is
  unpredictable.

  If the original pattern \code{x} has marks, then the replacement
  pattern \code{value} must also have marks, of the same type.

  Use the function \code{\link{unmark}} to remove marks from a
  marked point pattern.

  Use the function \code{\link{split.ppp}} to select those points
  in a marked point pattern which have a specified mark.
}
\seealso{
  \code{\link{subset.ppp}}.

  \code{\link{ppp.object}},
  \code{\link{owin.object}},
  \code{\link{unmark}},
  \code{\link{split.ppp}},
  \code{\link{cut.ppp}}
}
\section{Warnings}{
  The function does not check whether \code{i} is a subset of
  \code{x$window}. Nor does it check whether \code{value} lies
  inside \code{x$window}.
}
\examples{
 # Longleaf pines data
 \dontrun{
 plot(longleaf)
 }
 \testonly{
	longleaf <- longleaf[seq(1,longleaf$n,by=10)]
 }
 # adult trees defined to have diameter at least 30 cm
 adult <- (longleaf$marks >= 30)
 longadult <- longleaf[adult]
 \dontrun{
 plot(longadult)
 }
 # note that the marks are still retained.
 # Use unmark(longadult) to remove the marks
 
 # New Zealand trees data
 \dontrun{
 plot(nztrees)          # plot shows a line of trees at the far right
 abline(v=148, lty=2)   # cut along this line
 }
 nzw <- owin(c(0,148),c(0,95)) # the subwindow
 # trim dataset to this subwindow
 nzsub <- nztrees[nzw]
 \dontrun{
 plot(nzsub)
 }

 # Redwood data
 \dontrun{
 plot(redwood)
 }
 # Random thinning: delete 60\% of data
 retain <- (runif(redwood$n) < 0.4)
 thinred <- redwood[retain]
 \dontrun{
 plot(thinred)
 }
 # Scramble 60\% of data
 modif <- (runif(redwood$n) < 0.6)
 scramble <- function(x) { runifpoint(x$n, x$window) }
 redwood[modif] <- scramble(redwood[modif])

 # Lansing woods data - multitype points
 \testonly{
    lansing <- lansing[seq(1, lansing$n, length=100)]
 }

 # Hickory trees
  hicks <- split(lansing)$hickory

 # Trees in subwindow
  win <- owin(c(0.3, 0.6),c(0.2, 0.5))
  lsub <- lansing[win]

 # Scramble the locations of trees in subwindow, retaining their marks
  lansing[win] <- scramble(lsub) \%mark\% (lsub$marks)

 # Extract oaks only
 oaknames <- c("redoak", "whiteoak", "blackoak")
 oak <- lansing[marks(lansing) \%in\% oaknames, drop=TRUE]
 oak <- subset(lansing, marks \%in\% oaknames, drop=TRUE)

 # To clip or not to clip
 X <- runifpoint(25, letterR)
 B <- owin(c(2.2, 3.9), c(2, 3.5))
 opa <- par(mfrow=c(1,2))
 plot(X, main="X[B]")
 plot(X[B], border="red", cols="red", add=TRUE, show.all=TRUE, main="")
 plot(X, main="X[B, clip=TRUE]")
 plot(B, add=TRUE, lty=2)
 plot(X[B, clip=TRUE], border="blue", cols="blue", add=TRUE,
      show.all=TRUE, main="")
 par(opa)
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  and
  Rolf Turner \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
