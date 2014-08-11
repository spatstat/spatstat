\name{fasp.object}
\alias{fasp.object}
\title{Function Arrays for Spatial Patterns}
\description{
  A class \code{"fasp"} to represent a \dQuote{matrix}
  of functions, amenable to plotting as a matrix of plot panels.
}
\details{
  An object of this class is a convenient way of storing
  (and later plotting, editing, etc)
  a set of functions \eqn{f_{i,j}(r)}{f[i,j](r)} of a real argument \eqn{r},
  defined for each possible pair \eqn{(i,j)} of indices
  \eqn{1 \le i,j \le n}{1 <= i,j <= n}. We may think of this
  as a matrix or array of functions \eqn{f_{i,j}}{f[i,j]}.

  Function arrays are particularly useful in the 
  analysis of a multitype point pattern (a point pattern in which 
  the points are identified as belonging to separate types).
  We may want to compute a summary function for the points
  of type \eqn{i} only, for each of the possible types \eqn{i}.
  This produces a \eqn{1 \times m}{1 * m} array of functions.
  Alternatively we may compute a summary function
  for each possible pair of types \eqn{(i,j)}.
  This produces an \eqn{m \times m}{m *  m} array of functions.

  For multitype point patterns the command \code{\link{alltypes}}
  will compute arrays of summary functions for each possible
  type or for each possible pair of types.
  The function \code{\link{alltypes}} 
  returns an object of class \code{"fasp"}.

  An object of class \code{"fasp"} is a list containing at least the
  following components:

  \describe{
    \item{fns}{
      A list of data frames, each representing one of the functions.
    }
    \item{which}{
      A matrix representing the spatial arrangement of the
      functions. If \code{which[i,j] = k}
      then the function represented by \code{fns[[k]]} should be plotted
      in the panel at position \eqn{(i,j)}. If \code{which[i,j] = NA}
      then nothing is plotted in that position.
    }
    \item{titles}{
      A list of character strings, providing suitable plotting titles
      for the functions.
    }
    \item{default.formulae}{
      A list of default formulae for plotting each of the functions.
    }
    \item{title}{
      A character string, giving a default title for the array
      when it is plotted.
    }
  }
}
\section{Functions available}{
  There are methods for \code{plot}, \code{print} and \code{"["}
  for this class.

  The plot method displays the entire array of functions.
  The method \code{\link{[.fasp}} selects a sub-array using the natural
  indices \code{i,j}.

  The command \code{\link{eval.fasp}} can be used to apply
  a transformation to each function in the array,
  and to combine two arrays.
}
\seealso{
  \code{\link{alltypes}},
  \code{\link{plot.fasp}},
  \code{\link{[.fasp}},
  \code{\link{eval.fasp}}
}
\examples{
  # multitype point pattern
  data(amacrine)
  GG <- alltypes(amacrine, "G")
  plot(GG)

  # select the row corresponding to cells of type "on"
  Gon <- GG["on", ]
  plot(Gon)

  # extract the G function for i = "on", j = "off"
  Gonoff <- GG["on", "off", drop=TRUE]

  # Fisher variance stabilising transformation
  GGfish <- eval.fasp(asin(sqrt(GG)))
  plot(GGfish)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{attribute}
 
 
