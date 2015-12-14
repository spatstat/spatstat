\name{range.fv}
\alias{range.fv}
\alias{max.fv}
\alias{min.fv}
\title{
  Range of Function Values
}
\description{
  Compute the range, maximum, or minimum of the
  function values in a summary function.
}
\usage{
  \method{range}{fv}(\dots, na.rm = TRUE, finite = na.rm)

  \method{max}{fv}(\dots, na.rm = TRUE, finite = na.rm)

  \method{min}{fv}(\dots, na.rm = TRUE, finite = na.rm)
}
\arguments{
  \item{\dots}{
    One or more function value tables (objects of class \code{"fv"}
    representing summary functions) or other data.
  }
  \item{na.rm}{
    Logical. Whether to ignore \code{NA} values.
  }
  \item{finite}{
    Logical. Whether to ignore values that are
    infinite, \code{NaN} or \code{NA}.
  }
}
\details{
  These are methods for the generic \code{\link[base]{range}},
  \code{\link[base]{max}} and \code{\link[base]{min}}.
  They compute the range, maximum, and minimum of the \emph{function} values
  that would be plotted on the \eqn{y} axis by default.

  For more complicated calculations, use \code{\link{with.fv}}.
}
\value{
  Numeric vector of length 2.
}
\seealso{
  \code{\link{with.fv}}
}
\examples{
   G <- Gest(cells)
   range(G)
   max(G)
   min(G)
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  ,
  Rolf Turner \email{r.turner@auckland.ac.nz} and Ege Rubak \email{rubak@math.aau.dk}.
}
\keyword{spatial}
\keyword{math}
