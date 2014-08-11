\name{rshift}
\alias{rshift}
\title{Random Shift}
\description{
  Randomly shifts the points of a point pattern or
  line segment pattern. Generic.
}
\usage{
   rshift(X, \dots)
}
\arguments{
  \item{X}{Pattern to be subjected to a random shift.
    A point pattern (class \code{"ppp"}),
    a line segment pattern (class \code{"psp"})
    or an object of class \code{"splitppp"}.
  }
  \item{\dots}{
    Arguments controlling the generation of the
    random shift vector, or specifying which parts of the pattern
    will be shifted. 
  }
}
\value{
  An object of the same type as \code{X}.
}
\details{
  This operation applies a random shift (vector displacement) to
  the points in a point pattern,
  or to the segments in a line segment pattern.

  The argument \code{X} may be 
  \itemize{
    \item
    a point pattern
    (an object of class \code{"ppp"})
    \item
    a line segment pattern
    (an object of class \code{"psp"})
    \item
    an object of class \code{"splitppp"}
    (basically a list of point patterns, obtained from
    \code{\link{split.ppp}}).
  }
  The function \code{rshift} is generic, with
  methods for the three classes \code{"ppp"},
  \code{"psp"} and \code{"splitppp"}.

  See the help pages for these methods, \code{\link{rshift.ppp}},
  \code{\link{rshift.psp}} and
  \code{\link{rshift.splitppp}}, for further information.
}
\seealso{
  \code{\link{rshift.ppp}},
  \code{\link{rshift.psp}},
  \code{\link{rshift.splitppp}}
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}

