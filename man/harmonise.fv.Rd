\name{harmonise.fv}
\alias{harmonise.fv}
\alias{harmonize.fv}
\title{Make Function Tables Compatible}
\description{
  Convert several objects of class \code{"fv"} 
  to the same values of the function argument.
}
\usage{
\method{harmonise}{fv}(\dots, strict=FALSE)

\method{harmonize}{fv}(\dots, strict=FALSE)
}
\arguments{
  \item{\dots}{
    Any number of function tables (objects of class \code{"fv"}).
  }
  \item{strict}{
    Logical. If \code{TRUE}, a column of data will be deleted
    if columns of the same name do not appear in every object.
  }
}
\details{
  A function value table (object of class \code{"fv"}) is
  essentially a data frame giving the values of a function \eqn{f(x)}
  (or several alternative estimates of this value)
  at equally-spaced values of the function argument \eqn{x}.

  The command \code{\link{harmonise}} is generic. This is the
  method for objects of class \code{"fv"}.
  
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

  If \code{strict=TRUE}, each column of data will be retained only if
  a column of the same name appears in all of the arguments \code{\dots}.
  This ensures that the resulting objects are strictly compatible
  in the sense of \code{\link{compatible.fv}},
  and can be combined using \code{\link{eval.fv}}
  or \code{\link{collapse.fv}}.
  
  If \code{strict=FALSE} (the default), this does not occur,
  and then the resulting objects are \bold{not} guaranteed to be compatible
  in the sense of \code{\link{compatible.fv}}.
}
\value{
  A list, of length equal to the number of arguments \code{\dots},
  whose entries are objects of class \code{"fv"}.
  If the arguments were named (\code{name=value}) then the return value
  also carries these names.
}
\author{\adrian
  
  ,
  \rolf
  
  and \ege.
}
\examples{
   H <- harmonise(K=Kest(cells), G=Gest(cells))
   H
   \dontrun{
      ## generates a warning about duplicated columns
      try(cbind(H$K, H$G))
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
