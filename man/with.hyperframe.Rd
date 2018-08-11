\name{with.hyperframe}
\alias{with.hyperframe}
\title{Evaluate an Expression in Each Row of a Hyperframe}
\description{
  An expression, involving the names of columns in a hyperframe,
  is evaluated separately for each row of the hyperframe.
}
\usage{
\method{with}{hyperframe}(data, expr, ...,
                         simplify = TRUE,
                         ee = NULL, enclos=NULL)
}
\arguments{
  \item{data}{A hyperframe (object of class \code{"hyperframe"})
    containing data.
  }
  \item{expr}{An \R language expression to be evaluated.}
  \item{\dots}{Ignored.}
  \item{simplify}{
    Logical. If \code{TRUE}, the return value
    will be simplified to a vector whenever possible.
  }
  \item{ee}{
    Alternative form of \code{expr}, as an object of class
    \code{"expression"}.
  }
  \item{enclos}{
    An environment in which to search for objects that are
    not found in the hyperframe. Defaults to \code{\link{parent.frame}()}.
  }
}
\details{
  This function evaluates the expression \code{expr} in each row
  of the hyperframe \code{data}. It is a method for the generic
  function \code{\link{with}}.
  
  The argument \code{expr} should be an \R language expression
  in which each variable name is either the name of a column in the
  hyperframe \code{data}, or the name of an object in the parent frame
  (the environment in which \code{with} was called.)
  The argument \code{ee} can be used as an alternative
  to \code{expr} and should be an expression object (of
  class \code{"expression"}).
  
  For each row of \code{data}, the expression will be evaluated
  so that variables which are column names of \code{data} are
  interpreted as the entries for those columns in the current row.

  For example, if a hyperframe \code{h} has columns
  called \code{A} and \code{B}, then \code{with(h, A != B)} inspects
  each row of \code{data} in turn,
  tests whether the entries in columns \code{A} and \code{B} are
  equal, and returns the \eqn{n} logical values.
}
\value{
  Normally a list of length
  \eqn{n} (where \eqn{n} is the number of rows) containing the results
  of evaluating the expression for each row. 
  If \code{simplify=TRUE} and each result is a single atomic value,
  then the result is a vector or factor 
  containing the same values.
}
\author{\adrian
  
  
  and \rolf
  
}
\seealso{
  \code{\link{hyperframe}},
  \code{\link{plot.hyperframe}}
}
\examples{
  # generate Poisson point patterns with intensities 10 to 100
   H <- hyperframe(L=seq(10,100, by=10))
   X <- with(H, rpoispp(L))
}
\keyword{spatial}
\keyword{manip}
\keyword{programming}
