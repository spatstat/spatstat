\name{subset.hyperframe}
\alias{subset.hyperframe}
\title{
  Subset of Hyperframe Satisfying A Condition
}
\description{
  Given a hyperframe, return the subset specified by
  imposing a condition on each row, and optionally by choosing
  only some of the columns.
}
\usage{
\method{subset}{hyperframe}(x, subset, select, \dots)
}
\arguments{
  \item{x}{
    A hyperframe pattern (object of class \code{"hyperframe"}.
  }
  \item{subset}{
    Logical expression indicating which points are to be kept.
    The expression may involve the names of columns of \code{x}
    and will be evaluated by \code{\link{with.hyperframe}}.
  }
  \item{select}{
    Expression indicating which columns of marks should be kept.
  }
  \item{\dots}{
    Arguments passed to \code{\link{[.hyperframe}}
    such as \code{drop} and \code{strip}.
  }
}
\details{
  This is a method for the generic function \code{\link{subset}}.
  It extracts the subset of rows of \code{x}
  that satisfy the logical expression
  \code{subset}, and retains only the columns of \code{x} that are
  specified by the expression \code{select}. The result is always a
  hyperframe.

  The argument \code{subset} determines the subset of rows that
  will be extracted. It should be a logical expression.
  It may involve the names of columns of \code{x}.
  The default is to keep all points.

  The argument \code{select} determines which columns of \code{x}
  will be retained.
  It should be an expression involving the names of columns
  (which will be interpreted as integers representing the positions of
  these columns). For example if there are columns named
  \code{A} to \code{Z}, then \code{select=D:F} is a valid expression
  and means that columns \code{D}, \code{E} and \code{F} will be
  retained. Similarly \code{select=-(A:C)} is valid and means that columns
  \code{A} to \code{C} will be deleted.  
  The default is to retain all columns.

  Setting \code{subset=FALSE} will remove all the rows.
  Setting \code{select=FALSE} will remove all the columns.
  
  The result is always a hyperframe.
}
\value{
  A hyperframe.
}  
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/},
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\seealso{
  \code{\link[base]{subset}},
  \code{\link{[.hyperframe}}
}
\examples{
 a <- subset(flu, virustype=="wt")

 aa <- subset(flu, minnndist(pattern) > 10)

 aaa <- subset(flu, virustype=="wt", select = -pattern)
}
\keyword{spatial}
\keyword{manip}
