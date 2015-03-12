\name{by.ppp}
\alias{by.ppp}
\title{Apply a Function to a Point Pattern Broken Down by Factor}
\description{
  Splits a point pattern into sub-patterns, and applies the function to each
  sub-pattern.
}
\usage{
 \method{by}{ppp}(data, INDICES=marks(data), FUN, ...)
}
\arguments{
  \item{data}{Point pattern (object of class \code{"ppp"}).}
  \item{INDICES}{Grouping variable. Either a factor, a pixel image with
    factor values, or a tessellation.}
  \item{FUN}{Function to be applied to subsets of \code{data}.}
  \item{\dots}{Additional arguments to \code{FUN}.}
}
\details{
  This is a method for the generic function \code{\link{by}}
  for point patterns (class \code{"ppp"}).

  The point pattern \code{data} is first divided into subsets
  according to \code{INDICES}. Then the function \code{FUN}
  is applied to each subset. The results of each computation are
  returned in a list.

  The argument \code{INDICES} may be
  \itemize{
    \item
    a factor, of length equal to the number of points in \code{data}.
    The levels of \code{INDICES}
    determine the destination of each point in \code{data}.
    The \code{i}th point of \code{data} will be placed in the sub-pattern
    \code{split.ppp(data)$l} where \code{l = f[i]}.
    \item
    a pixel image (object of class \code{"im"}) with factor values.
    The pixel value of \code{INDICES}
    at each point of \code{data} will be used as the classifying variable.
    \item
    a tessellation (object of class \code{"tess"}).
    Each point of \code{data} will be classified according to
    the tile of the tessellation into which it falls.
  }
  If \code{INDICES} is missing, then \code{data} must be a multitype point pattern
  (a marked point pattern whose marks vector is a factor).
  Then the effect is that the points of each type
  are separated into different point patterns.
}
\value{
  A list (also of class \code{"anylist"} or \code{"solist"} as
  appropriate) containing the results returned
  from \code{FUN} for each of the subpatterns.
}
\seealso{
  \code{\link{ppp}},
  \code{\link{split.ppp}},
  \code{\link{cut.ppp}},
  \code{\link{tess}},
  \code{\link{im}}.
}
\examples{
  # multitype point pattern, broken down by type
  data(amacrine)
  by(amacrine, FUN=density)
  by(amacrine, FUN=function(x) { min(nndist(x)) } )

  # how to pass additional arguments to FUN
  by(amacrine, FUN=clarkevans, correction=c("Donnelly","cdf"))

  # point pattern broken down by tessellation
  data(swedishpines)
  tes <- quadrats(swedishpines, 5, 5)
  B <- by(swedishpines, tes, clarkevans, correction="Donnelly")
  unlist(lapply(B, as.numeric))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{methods}
\keyword{manip}
