\name{subset.psp}
\alias{subset.psp}
\title{
  Subset of Line Segment Satisfying A Condition
}
\description{
  Given a line segment pattern, return the subset
  of segments which satisfy a specified condition.
}
\usage{
\method{subset}{psp}(x, subset, select, drop=FALSE, \dots)
}
\arguments{
  \item{x}{
    A line segment pattern (object of class \code{"psp"}).
  }
  \item{subset}{
    Logical expression indicating which points are to be kept.
    The expression may involve the names of spatial coordinates
    of the segment endpoints
    (\code{x0}, \code{y0}, \code{x1}, \code{y1}), the \code{marks}, and
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
  It extracts the subset of \code{x} consisting of those segments
  that satisfy the logical expression
  \code{subset}, and retains only the columns of marks that are
  specified by the expression \code{select}. The result is
  always a line segment pattern, with the same window as \code{x}.

  The argument \code{subset} determines the subset that
  will be extracted. It should be a logical expression.
  It may involve the variable names
  \code{x0}, \code{y0}, \code{x1}, \code{y1}
  representing the Cartesian coordinates of the segment endpoints;
  the name \code{marks} representing the marks;
  and (if there is more than one column of marks)
  the names of individual columns of marks.
  The default is to keep all segments.

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

  The result is always a line segment pattern.
  To extract only some columns of marks as a data frame,
  use \code{subset(as.data.frame(x), ...)}
}
\value{
  A line segment pattern (object of class \code{"psp"})
  in the same spatial window as \code{x}. The result is a subset of \code{x},
  possibly with some columns of marks removed.
}  
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link[base]{subset}},
  \code{\link{[.psp}}.
}
\examples{
 plot(nbw.seg)
 plot(subset(nbw.seg, x0 < 500 & y0 < 800), add=TRUE, lwd=6)
 subset(nbw.seg, type == "island")
 subset(nbw.seg, type == "coast", select= -type)
 subset(nbw.seg, type \%in\% c("island", "coast"), select= FALSE)
}
\keyword{spatial}
\keyword{manip}
