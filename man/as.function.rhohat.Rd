\name{as.function.rhohat}
\alias{as.function.rhohat}
\title{
  Convert Function Table to Function
}
\description{
  Converts an object of class \code{"rhohat"} to an \R language function.
}
\usage{
  \method{as.function}{rhohat}(x, ..., value=".y", extrapolate=TRUE)
}
\arguments{
  \item{x}{
    Object of class \code{"rhohat"},
    produced by the function \code{\link{rhohat}}.
}
  \item{\dots}{
    Ignored.
}
  \item{value}{
    Optional. Character string or character vector selecting
    one or more of the columns of \code{x}
    for use as the function value. See Details.
  }
  \item{extrapolate}{
    Logical, indicating whether to extrapolate the function
    outside the domain of \code{x}. See Details.
  }
}
\details{
  An object of class \code{"rhohat"} is essentially a data frame
  of estimated values of the function \eqn{rho(x)}
  as described in the help file for \code{\link{rhohat}}. 
  
  Sometimes it is useful to convert the function value table
  to a \code{function} in the \R language. This is done by
  \code{as.function.rhohat}. It converts an object \code{x}
  of class \code{"rhohat"} to an \R function \code{f}.

  The command \code{as.function.rhohat} is a method for the generic command
  \code{\link{as.function}} for the class \code{"rhohat"}.
  
  If \code{f <- as.function(x)} then \code{f} is an \R function
  that accepts a numeric argument and returns a corresponding value
  for the summary function by linear interpolation between the values
  in the table \code{x}. 

  Argument values lying outside the range of the table
  yield an \code{NA} value (if \code{extrapolate=FALSE})
  or the function value at the nearest endpoint of the range
  (if \code{extrapolate = TRUE}).
  To apply different rules to the left and right extremes,
  use \code{extrapolate=c(TRUE,FALSE)} and so on.
  
  Typically the table \code{x} contains several columns of
  function values corresponding to different edge corrections.
  Auxiliary information for the table identifies one of these
  columns as the \emph{recommended value}. 
  By default, the values of the function \code{f <- as.function(x)}
  are taken from this column of recommended values. 
  This default can be changed using the argument \code{value},
  which can be a character string or character vector of names of
  columns of \code{x}. Alternatively \code{value} can be one of
  the abbreviations used by \code{\link{fvnames}}.

  If \code{value} specifies a single column of the table,
  then the result is a function \code{f(r)} with a single numeric
  argument \code{r} (with the same name as the orginal argument
  of the function table).

  If \code{value} specifies several columns of the table,
  then the result is a function \code{f(r,what)} 
  where \code{r} is the numeric argument and 
  \code{what} is a character string identifying the column of values
  to be used.

  The formal arguments of the resulting function
  are \code{f(r, what=value)}, which
  means that in a call to this function \code{f}, the permissible values
  of \code{what} are the entries of the original vector \code{value};
  the default value of \code{what} is the first entry of \code{value}.
  
}
\value{
  A \code{function(r)} or \code{function(r,what)}
  where \code{r} is the name of the original argument of the function table.
}
\author{
  \adrian
  and \rolf
}
\seealso{
   \code{\link{rhohat}}, 
   \code{\link{methods.rhohat}},
   \code{\link{as.function.fv}}. 
}
\examples{
  g <- rhohat(cells, "x")
  f <- as.function(g)
  f
  f(0.1)
}
\keyword{spatial}
\keyword{methods}

