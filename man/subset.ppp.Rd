\name{subset.ppp}
\alias{subset.ppp}
\alias{subset.pp3}
\alias{subset.lpp}
\alias{subset.ppx}
\title{
  Subset of Point Pattern Satisfying A Condition
}
\description{
  Given a point pattern, return the subset of points which satisfy
  a specified condition.
}
\usage{
\method{subset}{ppp}(x, subset, select, drop=FALSE, \dots)

\method{subset}{pp3}(x, subset, select, drop=FALSE, \dots)

\method{subset}{lpp}(x, subset, select, drop=FALSE, \dots)

\method{subset}{ppx}(x, subset, select, drop=FALSE, \dots)
}
\arguments{
  \item{x}{
    A point pattern (object of class \code{"ppp"},
    \code{"lpp"}, \code{"pp3"} or \code{"ppx"}).
  }
  \item{subset}{
    Logical expression indicating which points are to be kept.
    The expression may involve the names of spatial coordinates
    (\code{x}, \code{y}, etc), the \code{marks}, and
    (if there is more than one column of marks)
    the names of individual columns of marks.
    Missing values are taken as false. See Details.
  }
  \item{select}{
    Expression indicating which columns of marks should be kept.
    The \emph{names} of columns of marks can be used in this expression,
    and will be treated as if they were column indices.
    See Details.
  }
  \item{drop}{
    Logical value indicating whether to remove unused levels
    of the marks, if the marks are a factor.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  This is a method for the generic function \code{\link{subset}}.
  It extracts the subset of points of \code{x}
  that satisfy the logical expression
  \code{subset}, and retains only the columns of marks that are
  specified by the expression \code{select}. The result is always a point
  pattern, with the same window as \code{x}.

  The argument \code{subset} determines the subset of points that
  will be extracted. It should be a logical expression.
  It may involve the variable names
  \code{x} and \code{y} representing the Cartesian coordinates;
  the names of other spatial coordinates or local coordinates;
  the name \code{marks} representing the marks;
  and (if there is more than one column of marks)
  the names of individual columns of marks.
  The default is to keep all points.

  The argument \code{select} determines which columns of marks
  will be retained (if there are several columns of marks).
  It should be an expression involving the names of columns of marks
  (which will be interpreted as integers representing the positions of
  these columns). For example if there are columns of marks named
  \code{A} to \code{Z}, then \code{select=D:F} is a valid expression
  and means that columns \code{D}, \code{E} and \code{F} will be
  retained. Similarly \code{select=-(A:C)} is valid and means that columns
  \code{A} to \code{C} will be deleted.  
  The default is to retain all columns.

  Setting \code{subset=FALSE} will produce an empty point pattern
  (i.e. containing zero points) in the same window as \code{x}.
  Setting \code{select=FALSE} or \code{select= -marks} will
  remove all the marks from \code{x}.

  The argument \code{drop} determines whether to remove
  unused levels of a factor, if the resulting point pattern is multitype
  (i.e. the marks are a factor) or if the marks are a data frame
  in which some of the columns are factors.

  The result is always a point pattern, of the same class as \code{x}.
  Spatial coordinates (and local
  coordinates) are always retained. To extract only some
  columns of marks or coordinates as a data frame,
  use \code{subset(as.data.frame(x), ...)}
}
\section{Other kinds of subset arguments}{
  Alternatively the argument \code{subset} can be any kind of subset index
  acceptable to \code{\link{[.ppp}}, \code{\link{[.pp3}},
  \code{\link{[.ppx}}. This argument selects which points of \code{x}
  will be retained.

  \bold{Warning:} if the argument \code{subset} is
  a window, this is interpreted as specifying the subset of points that fall
  inside that window, but the resulting point pattern has the same window
  as the original pattern \code{x}.
}
\value{
  A point pattern of the same class as \code{x}, in the same
  spatial window as \code{x}. The result is a subset of \code{x},
  possibly with some columns of marks removed.
}  
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link[base]{subset}},
  
  \code{\link{[.ppp}},
  \code{\link{[.pp3}},
  \code{\link{[.lpp}},
  \code{\link{[.ppx}}
}
\examples{
 plot(subset(cells, x > 0.5))

 subset(amacrine, marks == "on")

 subset(amacrine, marks == "on", drop=TRUE)

 subset(redwood, nndist(redwood) > 0.04)

 subset(finpines, select=height)

 subset(finpines, diameter > 2, height)

 subset(nbfires, year==1999 & ign.src == "campfire",
                 select=cause:fnl.size)

 v <- subset(chicago, x + y > 1100 & marks == "assault")

 vv <- subset(chicago, x + y > 1100 & marks == "assault", drop=TRUE)

 a <- subset(rpoispp3(40), z > 0.5)
}
\keyword{spatial}
\keyword{manip}
