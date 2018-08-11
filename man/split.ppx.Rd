\name{split.ppx}
\alias{split.ppx}
\title{Divide Multidimensional Point Pattern into Sub-patterns}
\description{
  Divides a multidimensional point pattern into several sub-patterns,
  according to their marks, or according to any user-specified grouping.
}
\usage{
  \method{split}{ppx}(x, f = marks(x), drop=FALSE, un=NULL, \dots)
}
\arguments{
  \item{x}{
    A multi-dimensional point pattern.
    An object of class \code{"ppx"}.
  }
  \item{f}{
    Data determining the grouping. Either a factor, 
    or the name of one of the columns of marks.
  }
  \item{drop}{
    Logical. Determines whether empty groups will
    be deleted.
  }
  \item{un}{
    Logical. Determines whether the resulting subpatterns will be unmarked
    (i.e. whether marks will be removed	from the points in each
    subpattern). 
  }
  \item{\dots}{
    Other arguments are ignored.
  }
} 
\value{
  A list of point patterns.
  The components of the list are named by the levels of \code{f}.
  The list also has the class \code{"splitppx"} and \code{"anylist"}.
}
\details{
  The generic command \code{\link[base]{split}} allows a dataset to be separated
  into subsets according to the value of a grouping variable.
  
  The function \code{split.ppx} is a method for the generic
  \code{\link[base]{split}} for the class \code{"ppx"} of multidimensional
  point patterns. It divides up the points of the point pattern \code{x}
  into several sub-patterns according to the values of \code{f}.
  The result is a list of point patterns.
  
  The argument \code{f} may be
  \itemize{
    \item
    a factor, of length equal to the number of points in \code{x}.
    The levels of \code{f}
    determine the destination of each point in \code{x}.
    The \code{i}th point of \code{x} will be placed in the sub-pattern
    \code{split.ppx(x)$l} where \code{l = f[i]}.
    \item
    a character string, matching the name of one of the columns of
    marks, if \code{marks(x)} is a data frame. This column should
    be a factor.
  }
  If \code{f} is missing, then it will be determined by the
  marks of the point pattern. The pattern \code{x} can be either
  \itemize{
    \item 
    a multitype point pattern
    (a marked point pattern whose marks vector is a factor).
    Then \code{f} is taken to be the marks vector.
    The effect is that the points of each type
    are separated into different point patterns.
    \item 
    a marked point pattern with a data frame or hyperframe
    of marks, containing at least one
    column that is a factor. The first such column will be used to
    determine the splitting factor \code{f}.
  }

  Some of the sub-patterns created by the split
  may be empty. If \code{drop=TRUE}, then empty sub-patterns will
  be deleted from the list. If \code{drop=FALSE} then they are retained.

  The argument \code{un} determines how to handle marks 
  in the case where \code{x} is a marked point pattern.
  If \code{un=TRUE} then the marks of the 
  points will be discarded when they are split into groups,
  while if \code{un=FALSE} then the marks will be retained.

  If \code{f} and \code{un} are both missing,
  then the default is \code{un=TRUE} for multitype point patterns
  and \code{un=FALSE} for marked point patterns with a data frame of
  marks.
  
  The result of \code{split.ppx} has class \code{"splitppx"}
  and \code{"anylist"}. There are methods for \code{print},
  \code{summary} and \code{plot}.
}
\seealso{
  \code{\link{ppx}},
  \code{\link{plot.anylist}}
}
\examples{
   df <- data.frame(x=runif(4),y=runif(4),t=runif(4),
                    age=rep(c("old", "new"), 2),
                    size=runif(4))
   X <- ppx(data=df, coord.type=c("s","s","t","m","m"))
   X
   split(X)
}

\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{methods}
\keyword{manip}
