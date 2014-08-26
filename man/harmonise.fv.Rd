\name{harmonise.fv}
\alias{harmonise.fv}
\alias{harmonize.fv}
\title{Make Function Tables Compatible}
\description{
  Convert several objects of class \code{"fv"} 
  to the same values of the function argument.
}
\usage{
harmonise.fv(\dots)
harmonize.fv(\dots)
}
\arguments{
  \item{\dots}{
    Any number of function tables (objects of class \code{"fv"}).
  } 
}
\details{
  A function value table (object of class \code{"fv"}) is
  essentially a data frame giving the values of a function \eqn{f(x)}
  (or several alternative estimates of this value)
  at equally-spaced values of the function argument \eqn{x}.

  This command makes any number of \code{"fv"} objects compatible,
  in the loose sense that they have the same sequence of values of
  \eqn{x}. They can then be combined by \code{\link{cbind.fv}},
  but not necessarily by \code{\link{eval.fv}}.

  All arguments \code{\dots} must be function value tables
  (objects of class \code{"fv"}). 
  The result will be a list, of length equal to the number of
  arguments \code{\dots}, containing new versions of each of these functions,
  converted to a common sequence of \eqn{x} values.
  If the arguments were named (\code{name=value}) then the return value
  also carries these names.

  The range of \eqn{x} values in the resulting functions
  will be the intersection of the ranges of \eqn{x} values
  in the original functions. 
  The spacing of \eqn{x} values in the resulting functions
  will be the finest (narrowest) of the spacings of the
  \eqn{x} values in the original functions.
  Function values are interpolated using \code{\link[stats]{approxfun}}.

  The resulting objects are \bold{not} guaranteed to be compatible
  in the strict sense of \code{\link{compatible.fv}},
  so they cannot necessarily be combined using \code{\link{eval.fv}}
  or \code{\link{collapse.fv}}.
}
\value{
  A list, of length equal to the number of arguments \code{\dots},
  whose entries are objects of class \code{"fv"}.
  If the arguments were named (\code{name=value}) then the return value
  also carries these names.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/},
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak.
}
\examples{
   H <- harmonise.fv(K=Kest(cells), G=Gest(cells))
   H
   \dontrun{
      ## generates a warning about duplicated columns
      cbind(H$K, H$G)
   }
}
\seealso{
  \code{\link{fv.object}},
  \code{\link{cbind.fv}},
  \code{\link{eval.fv}},
  \code{\link{compatible.fv}}
}
\keyword{spatial}
\keyword{manip}
