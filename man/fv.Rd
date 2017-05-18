\name{fv}
\alias{fv}
\title{
  Create a Function Value Table
}
\description{
  Advanced Use Only.
  This low-level function creates an object of class \code{"fv"}
  from raw numerical data. 
}
\usage{
fv(x, argu = "r", ylab = NULL, valu, fmla = NULL, alim = NULL,
   labl = names(x), desc = NULL, unitname = NULL, fname = NULL, yexp = ylab)
}
\arguments{
  \item{x}{
    A data frame with at least 2 columns containing the
    values of the function argument and the corresponding values
    of (one or more versions of) the function.
  }
  \item{argu}{
    String. The name of the column of \code{x} that contains
    the values of the function argument.
  }
  \item{ylab}{
    Either \code{NULL}, or an \R language expression
    representing the mathematical name of the
    function. See Details.
  }
  \item{valu}{
    String. The name of the column of \code{x} that should be taken
    as containing the function values, in cases where a single column
    is required.
  }
  \item{fmla}{
    Either \code{NULL}, or a \code{formula} specifying the default
    plotting behaviour. See Details.
  }
  \item{alim}{
    Optional.
    The default range of values of the function argument for which the
    function will be plotted.
    Numeric vector of length 2.
  }
  \item{labl}{
    Optional. Plot labels for the columns of \code{x}.
    A vector of strings, with one entry for each column of \code{x}.
  }
  \item{desc}{
    Optional. Descriptions of the columns of \code{x}.
    A vector of strings, with one entry for each column of \code{x}.
  }
  \item{unitname}{
    Optional. Name of the unit (usually a unit of length)
    in which the function argument is
    expressed. Either a single character string,
    or a vector of two character strings giving the
    singular and plural forms, respectively.
  }
  \item{fname}{
    Optional.
    The name of the function itself. A character string.
  }
  \item{yexp}{
    Optional. Alternative form of \code{ylab}
    more suitable for annotating an axis of the plot.
    See Details.
  }
}
\details{
  This documentation is provided
  for experienced programmers who want to modify the internal
  behaviour of \pkg{spatstat}. Other users please see \code{\link{fv.object}}.

  The low-level function \code{fv} is used to create an object of
  class \code{"fv"} from raw numerical data.

  The data frame \code{x} contains the numerical data.
  It should have one column
  (typically but not necessarily named \code{"r"})
  giving the values of the function argument for which
  the function has been evaluated; and at least one other column,
  containing the corresponding values of the function.

  Typically there is more than one column of function values.
  These columns typically give the values of different versions or estimates
  of the same function,
  for example, different estimates of the \eqn{K} function
  obtained using different edge corrections.
  However they may also contain the values of related functions
  such as the derivative or hazard rate.

  \code{argu} specifies the name of the column of
  \code{x} that contains the values of the function argument
  (typically \code{argu="r"} but this is not compulsory).

  \code{valu} specifies the name of another column
  that contains the \sQuote{recommended} estimate of the function.
  It will be used to provide function values in those situations where
  a single column of data is required. For example,
  \code{\link{envelope}} computes its simulation envelopes
  using the recommended value of the summary function.

  \code{fmla} specifies the default plotting behaviour.
  It should be a formula, or a string that can be converted to a
  formula. Variables in the formula are names of columns of \code{x}.
  See \code{\link{plot.fv}} for the interpretation of this
  formula.

  \code{alim} specifies the recommended range of the
  function argument. This is used in situations where statistical
  theory or statistical practice indicates that the computed
  estimates of the function are not trustworthy outside a certain
  range of values of the function argument. By default,
  \code{\link{plot.fv}} will restrict the plot to this range.

  \code{fname} is a string giving the name of the function itself.
  For example, the \eqn{K} function would have \code{fname="K"}.
  
  \code{ylab} is a mathematical expression
  for the function value, used when labelling an axis
  of the plot, or when printing a description of the
  function. It should be an \R language object. 
  For example the \eqn{K} function's mathematical name \eqn{K(r)} is rendered
  by \code{ylab=quote(K(r))}. 

  If \code{yexp} is present, then \code{ylab} will be
  used only for printing, and \code{yexp} will be used for
  annotating axes in a plot. (Otherwise \code{yexp} defaults to \code{ylab}).
  For example the cross-type \eqn{K} function
  \eqn{K_{1,2}(r)}{K[1,2](r)} is rendered by something like
  \code{ylab=quote(Kcross[1,2](r))}
  and 
  \code{yexp=quote(Kcross[list(1,2)](r))}
  to get the most satisfactory behaviour.

  (A useful tip: use \code{\link{substitute}} instead of
  \code{\link{quote}} to insert values of variables into an expression,
  e.g. \code{substitute(Kcross[i,j](r), list(i=42,j=97))}
  yields the same as \code{quote(Kcross[42, 97](r))}.)

  \code{labl} is a character vector specifying plot labels
  for each column of \code{x}. These labels will appear on the
  plot axes (in non-default plots), legends and printed output.
  Entries in \code{labl}
  may contain the string \code{"\%s"} which will be replaced
  by \code{fname}. For example the border-corrected estimate
  of the \eqn{K} function has label \code{"\%s[bord](r)"} which
  becomes \code{"K[bord](r)"}.

  \code{desc} is a character vector containing intelligible
  explanations of each column of \code{x}. Entries in
  \code{desc} may contain the string \code{"\%s"} which will be replaced
  by \code{ylab}. For example the border correction estimate of the
  \eqn{K} function has description \code{"border correction estimate of \%s"}.

}
\value{
  An object of class \code{"fv"}, see \code{\link{fv.object}}.
}
\seealso{
  See \code{\link{plot.fv}} for plotting an \code{"fv"} object.

  See \code{\link{as.function.fv}} to convert an \code{"fv"} object
  to an \R function.
  
  Use \code{\link{cbind.fv}} to combine several \code{"fv"} objects.
  Use \code{\link{bind.fv}} to glue additional columns onto an existing
  \code{"fv"} object.

  Use \code{\link{range.fv}} to compute the range of \eqn{y} values
  for a function, and \code{\link{with.fv}} for more complicated
  calculations.
  
  The functions \code{fvnames}, \code{fvnames<-} allow the user to
  use standard abbreviations to refer to columns of an \code{"fv"} object.
  
  \emph{Undocumented} functions for modifying an \code{"fv"} object
  include \code{tweak.fv.entry} and \code{rebadge.fv}.
}
\author{\adrian
  
  
  and \rolf
  
}
\examples{
  df <- data.frame(r=seq(0,5,by=0.1))
  df <- transform(df, a=pi*r^2, b=3*r^2)
  X <- fv(df, "r", quote(A(r)),
              "a", cbind(a, b) ~ r,
              alim=c(0,4),
              labl=c("r", "\%s[true](r)", "\%s[approx](r)"),
              desc=c("radius of circle",
                     "true area \%s",
                     "rough area \%s"),
              fname="A")
  X
}
\keyword{spatial}
\keyword{classes}

