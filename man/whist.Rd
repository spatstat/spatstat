\name{whist}
\alias{whist}
\title{
  Weighted Histogram
}
\description{
  Computes the weighted histogram of a set of observations
  with a given set of weights.
}
\usage{
whist(x, breaks, weights = NULL)
}
\arguments{
  \item{x}{
    Numeric vector of observed values.
  }
  \item{breaks}{
    Vector of breakpoints for the histogram. 
  }
  \item{weights}{
    Numeric vector of weights for the observed values. 
  }
}
\details{
  This low-level function computes (but does not plot) the weighted
  histogram of a vector of observations \code{x} using a given
  vector of \code{weights}.

  The arguments \code{x} and \code{weights} should be numeric vectors of
  equal length. They may include \code{NA} or infinite values.

  The argument \code{breaks} should be a numeric vector whose entries
  are strictly increasing. These values define the boundaries between the
  successive histogram cells.
  The breaks \emph{do not} have to span the range
  of the observations.

  There are \code{N-1} histogram cells, where \code{N = length(breaks)}.
  An observation \code{x[i]} falls in the \code{j}th cell if
  \code{breaks[j] <= x[i] < breaks[j+1]} (for \code{j < N-1})
  or
  \code{breaks[j] <= x[i] <= breaks[j+1]} (for \code{j = N-1}).
  The weighted histogram value \code{h[j]} for the \code{j}th cell is
  the sum of \code{weights[i]} for all observations \code{x[i]} that
  fall in the cell.
  
  Note that, in contrast to the function \code{\link{hist}},
  the function \code{whist} does not require the breakpoints to span the
  range of the observations \code{x}. Values of \code{x} that fall outside the
  range of \code{breaks} are handled separately; their total weight
  is returned as an attribute of the histogram. 
}
\value{
  A numeric vector of length \code{N-1} containing the
  histogram values, where \code{N = length(breaks)}.

  The return value also has attributes \code{"low"} and \code{"high"}
  giving the total weight of all observations that are less than
  the lowest breakpoint, or greater than the highest breakpoint,
  respectively.
}
\examples{
  x <- rnorm(100)
  b <- seq(-1,1,length=21)
  w <- runif(100)
  whist(x,b,w)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
  with thanks to Peter Dalgaard.
}
\keyword{arith}
