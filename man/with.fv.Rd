\name{with.fv}
\alias{with.fv}
\title{Evaluate an Expression in a Function Table}
\description{
  Evaluate an R expression in a 
  function value table (object of class \code{"fv"}).
}
\usage{
\method{with}{fv}(data, expr, ..., drop = TRUE, enclos=NULL)
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
  \item{drop}{Logical value. If the result of evaluating the expression
    \code{expr} is a vector (rather than a matrix or data frame)
    then the result will be returned as a vector if \code{drop=TRUE}.
    Otherwise it will be returned as another function value table
    (object of class \code{"fv"}).
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

  The expression should be capable of handling vectors and matrices. 
  If the result of evaluating the expression is a matrix or data frame,
  then it is returned as a new function value table
  (object of class \code{"fv"}). If the result of evaluation is
  a vector and \code{drop=TRUE} then the result is returned as a vector.

  To perform calculations involving \emph{several} objects of
  class \code{"fv"}, use \code{\link{eval.fv}}.
}
\value{
  Either a function value table (object of class \code{"fv"})
  or a vector.
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
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
\keyword{programming}
