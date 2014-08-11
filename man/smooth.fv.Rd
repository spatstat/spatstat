\name{smooth.fv}
\alias{smooth.fv}
\alias{Smooth.fv}
\title{
  Apply Smoothing to Function Values
}
\description{
  Applies smoothing to the values
  in selected columns of a function value table.
}
\usage{
smooth.fv(x, which = "*", ...,
          method=c("smooth.spline", "loess"),
          xinterval=NULL)

\method{Smooth}{fv}(X, which = "*", ...,
          method=c("smooth.spline", "loess"),
          xinterval=NULL)
}
\arguments{
  \item{x,X}{
    Values to be smoothed.
    A function value table (object of class \code{"fv"},
    see \code{\link{fv.object}}).
  }
  \item{which}{
    Character vector identifying which columns of the table
    should be smoothed. Either a vector containing names
    of columns, or one of the wildcard strings \code{"*"} or \code{"."}
    explained below.
  }
  \item{\dots}{
    Extra arguments passed to \code{\link[stats]{smooth.spline}}
    or \code{\link[stats]{loess}} to control the smoothing.
  }
  \item{method}{
    Smoothing algorithm. A character string, partially matched
    to either \code{"smooth.spline"} or \code{"loess"}.
  }
  \item{xinterval}{
    Optional. Numeric vector of length 2 specifying a range of
    \eqn{x} values. Smoothing will be performed only on the part of the
    function corresponding to this range.
  }
}
\details{
  The command \code{smooth.fv} or \code{Smooth.fv}
  applies smoothing to the function values in
  a function value table (object of class \code{"fv"}).

  \code{Smooth.fv} is a method for the generic function
  \code{\link{Smooth}}.
  
  The smoothing is performed either by 
  \code{\link[stats]{smooth.spline}} or by
  \code{\link[stats]{loess}}.

  Smoothing is applied to every column
  (or to each of the selected columns) of function values in turn,
  using the function argument as the \eqn{x} coordinate
  and the selected column as the \eqn{y} coordinate.
  The original function values are then replaced by the corresponding
  smooth interpolated function values.

  The optional argument \code{which} specifies which of the
  columns of function values in \code{x} will be smoothed.
  The default (indicated by the wildcard \code{which="*"})
  is to smooth all function values, i.e.\ all columns except the
  function argument. Alternatively \code{which="."} designates
  the subset of function values that are displayed in the default plot.
  Alternatively \code{which} can be a character vector containing the
  names of columns of \code{x}.

  If the argument \code{xinterval} is given, then
  smoothing will be performed only in the specified range 
  of \eqn{x} values. 
}
\value{
  Another function value table (object of class \code{"fv"})
  of the same format.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{Smooth}},
  \code{\link{with.fv}},
  \code{\link{fv.object}},
  \code{\link[stats]{smooth.spline}},
  \code{\link[stats]{smooth.spline}}
}
\examples{
   data(cells)
   G <- Gest(cells)
   plot(G)
   plot(Smooth(G, df=9), add=TRUE)
}
\keyword{spatial}
\keyword{nonparametric}

