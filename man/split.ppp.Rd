\name{split.ppp}
\alias{split.ppp}
\alias{split<-.ppp}
\title{Divide Point Pattern into Sub-patterns}
\description{
  Divides a point pattern into several sub-patterns,
  according to their marks, or according to any user-specified grouping.
}
\usage{
  \method{split}{ppp}(x, f = marks(x), drop=FALSE, un=NULL, reduce=FALSE, \dots)
  \method{split}{ppp}(x, f = marks(x), drop=FALSE, un=missing(f), \dots) <- value
}
\arguments{
  \item{x}{
    A two-dimensional point pattern.
    An object of class \code{"ppp"}.
  }
  \item{f}{
    Data determining the grouping. Either a factor,
    a logical vector,
    a pixel image with factor values, a tessellation, a window,
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
  \item{reduce}{
    Logical. Determines whether to delete the column of marks used to split the
    pattern, when the marks are a data frame.
  }
  \item{\dots}{
    Other arguments are ignored.
  }
  \item{value}{
    List of point patterns.
  }
} 
\value{
  The value of \code{split.ppp} is a list of point patterns.
  The components of the list are named by the levels of \code{f}.
  The list also has the class \code{"splitppp"}.

  The assignment form \code{split<-.ppp} returns the updated
  point pattern \code{x}. 
}
\details{
  The function \code{split.ppp}
  divides up the points of the point pattern \code{x}
  into several sub-patterns according to the values of \code{f}.
  The result is a list of point patterns.
  
  The argument \code{f} may be
  \itemize{
    \item
    a factor, of length equal to the number of points in \code{x}.
    The levels of \code{f}
    determine the destination of each point in \code{x}.
    The \code{i}th point of \code{x} will be placed in the sub-pattern
    \code{split.ppp(x)$l} where \code{l = f[i]}.
    \item
    a pixel image (object of class \code{"im"}) with factor values.
    The pixel value of \code{f}
    at each point of \code{x} will be used as the classifying variable.
    \item
    a tessellation (object of class \code{"tess"}).
    Each point of \code{x} will be classified according to
    the tile of the tessellation into which it falls.
    \item
    a window (object of class \code{"owin"}).
    Each point of \code{x} will be classified according to
    whether it falls inside or outside this window.
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
    a marked point pattern with a data frame of marks, containing at least one
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

  If the marks of \code{x} are a data frame, then 
  \code{split(x, reduce=TRUE)} will discard only the column of marks
  that was used to split the pattern. This applies only when
  the argument \code{f} is missing.

  The result of \code{split.ppp} has class \code{"splitppp"}
  and can be plotted using \code{\link{plot.splitppp}}.

  The assignment function \code{split<-.ppp} 
  updates the point pattern \code{x} so that
  it satisfies \code{split(x, f, drop, un) = value}. The argument \code{value}
  is expected to be a list of point patterns, one for each level of
  \code{f}. These point patterns are expected to be compatible with the
  type of data in the original pattern \code{x}.

  Splitting can also be undone by the function
  \code{\link{superimpose}},
  but this typically changes the ordering of the data.
}
\seealso{
  \code{\link{cut.ppp}},
  \code{\link{plot.splitppp}},
  \code{\link{superimpose}},
  \code{\link{im}},
  \code{\link{tess}},
  \code{\link{ppp.object}}
}
\examples{

# (1) Splitting by marks

# Multitype point pattern: separate into types
 data(amacrine)
 u <- split(amacrine)

# plot them
 plot(split(amacrine))

# the following are equivalent:
 amon <- split(amacrine)$on
 amon <- unmark(amacrine[amacrine$marks == "on"])
   
# the following are equivalent:
 amon <- split(amacrine, un=FALSE)$on
 amon <- amacrine[amacrine$marks == "on"]
   
# Scramble the locations of the 'on' cells
 u <- split(amacrine)
 u$on <- runifpoint(amon$n, amon$window)
 split(amacrine) <- u

# Point pattern with continuous marks
 data(longleaf)
 \testonly{
	# smaller dataset
	longleaf <- longleaf[seq(1, longleaf$n, by=80)]
 }
 # cut the range of tree diameters into three intervals
 # using cut.ppp
 long3 <- cut(longleaf, breaks=3)
 # now split them
 long3split <- split(long3)

# (2) Splitting by a factor

# Unmarked point pattern
  data(swedishpines)
# cut & split according to nearest neighbour distance
  f <- cut(nndist(swedishpines), 3)
  u <- split(swedishpines, f)

# (3) Splitting over a tessellation
   tes <- tess(xgrid=seq(0,96,length=5),ygrid=seq(0,100,length=5))
   v <- split(swedishpines, tes)


# (4) how to apply an operation to selected points:
#  split into components, transform desired component, then un-split
#  e.g. apply random jitter to 'on' points only
  X <- amacrine
  Y <- split(X)
  Y$on <- rjitter(Y$on, 0.1)
  split(X) <- Y
}

\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{methods}
\keyword{manip}
