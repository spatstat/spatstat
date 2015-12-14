\name{with.fv}
\alias{with.fv}
\title{Evaluate an Expression in a Function Table}
\description{
  Evaluate an R expression in a 
  function value table (object of class \code{"fv"}).
}
\usage{
\method{with}{fv}(data, expr, ..., fun = NULL, enclos=NULL)
}
\arguments{
  \item{data}{A function value table (object of class \code{"fv"})
    in which the expression will be evaluated.
  }
  \item{expr}{The expression to be evaluated. An \R language
    expression, which may involve the names of columns in \code{data},
    the special abbreviations \code{.}, \code{.x} and \code{.y},
    and global constants or functions.
  }
  \item{\dots}{Ignored.}
  \item{fun}{Logical value, specifying whether the result
    should be interpreted as another function (\code{fun=TRUE})
    or simply returned as a numeric vector or array (\code{fun=FALSE}).
    See Details.
  }
  \item{enclos}{
    An environment in which to search for variables that are
    not found in \code{data}. Defaults to \code{\link{parent.frame}()}.
  }
}
\details{
  This is a method for the generic command \code{\link{with}}
  for an object of class \code{"fv"} (function value table).

  An object of class \code{"fv"} is a convenient way of storing and
  plotting several different estimates of the same function. It is
  effectively a data frame with extra attributes. 
  See \code{\link{fv.object}} for further explanation.

  This command makes it possible to perform computations that involve
  different estimates of the same function. For example we use it to compute
  the arithmetic difference between two different edge-corrected
  estimates of the \eqn{K} function of a point pattern.

  The argument \code{expr} should be an \R language expression. The expression
  may involve
  \itemize{
    \item the name of any column in \code{data}, referring to
    one of the estimates of the function;
    \item the symbol \code{.} which stands for all
    the available estimates of the function;
    \item the symbol \code{.y} which stands for the recommended
    estimate of the function (in an \code{"fv"} object, one of the
    estimates is always identified as the recommended estimate);
    \item the symbol \code{.x} which stands for the argument of the function;
    \item global constants or functions.
  }
  See the Examples.
  The expression should be capable of handling
  vectors and matrices.

  The interpretation of the argument \code{fun} is as follows:
  \itemize{
    \item
    If \code{fun=FALSE}, the result of evaluating the expression
    \code{expr} will be returned as a numeric vector, matrix or 
    data frame.
    \item
    If \code{fun=TRUE}, then the result of evaluating \code{expr}
    will be interpreted as containing the values of a new function.
    The return value will be an object of class \code{"fv"}.
    (This can only happen if the result has the right dimensions.)
    \item 
    The default is \code{fun=TRUE} if the result of evaluating
    \code{expr} has more than one column, and \code{fun=FALSE} otherwise.
  }
  
  To perform calculations involving \emph{several} objects of
  class \code{"fv"}, use \code{\link{eval.fv}}.
}
\value{
  A function value table (object of class \code{"fv"})
  or a numeric vector or data frame.
}
\seealso{
  \code{\link{with}},
  \code{\link{fv.object}},
  \code{\link{eval.fv}},
  \code{\link{Kest}}
}
\examples{
  # compute 4 estimates of the K function
  X <- rpoispp(42)
  K <- Kest(X)
  plot(K)

  # derive 4 estimates of the L function L(r) = sqrt(K(r)/pi)
  L <- with(K, sqrt(./pi))
  plot(L)

  # compute 4 estimates of V(r) = L(r)/r
  V <- with(L, ./.x)
  plot(V)

  # compute the maximum absolute difference between
  # the isotropic and translation correction estimates of K(r)
  D <- with(K, max(abs(iso - trans)))
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{manip}
\keyword{programming}
