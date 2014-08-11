\name{smooth.fv}
\alias{smooth.fv}
\title{
  Apply Smoothing to Function Values
}
\description{
  Applies spline smoothing to the values
  in selected columns of a function value table.
}
\usage{
smooth.fv(x, which = "*", ...)
}
\arguments{
  \item{x}{
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
    Extra arguments passed to \code{\link{smooth.spline}} to
    control the smoothing.
  }
}
\details{
  \code{\link{smooth.spline}} is applied 
  to each of the selected columns in turn
  (using the function argument as the \eqn{x} coordinate
  and the selected column as the \eqn{y} coordinate).
  The original function values are then replaced by the corresponding
  smooth interpolated function values.

  The argument \code{which} specifies which of the
  columns of function values in \code{x} will be smoothed.
  The default (indicated by the wildcard \code{which="*"})
  is to smooth all function values, i.e.\ all columns except the
  function argument. Alternatively \code{which="."} designates
  the subset of function values that are displayed in the default plot.
  Alternatively \code{which} can be a character vector containing the
  names of columns of \code{x}.
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
  \code{\link{with.fv}},
  \code{\link{fv.object}},
  \code{\link{smooth.spline}}
}
\examples{
   data(cells)
   G <- Gest(cells)
   plot(G)
   plot(smooth.fv(G, df=9), add=TRUE)
}
\keyword{spatial}
\keyword{nonparametric}

