\name{summary.splitppp}
\alias{summary.splitppp}
\title{Summary of a Split Point Pattern}
\description{
  Prints a useful summary of a split point pattern.
}
\usage{
 \method{summary}{splitppp}(object, \dots)
}
\arguments{
  \item{object}{
    Split point pattern (object of class \code{"splitppp"}, effectively a
    list of point patterns, usually created by \code{\link{split.ppp}}).
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  This is a method for the generic function \code{\link{summary}}.

  An object of the class \code{"splitppp"} is effectively a list
  of point patterns (objects of class \code{"ppp"}) representing
  different sub-patterns of an original point pattern.

  This function extracts a useful summary of each of the sub-patterns.
}
\seealso{
  \code{\link{summary}},
  \code{\link{split}},
  \code{\link{split.ppp}}
}
\examples{
  data(amacrine)      # multitype point pattern
  summary(split(amacrine))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{methods}

