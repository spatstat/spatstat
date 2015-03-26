\name{relevel.im}
\alias{relevel.im}
\alias{relevel.ppp}
\alias{relevel.ppx}
\title{
  Reorder Levels of a Factor-Valued Image or Pattern
}
\description{
  For a pixel image with factor values, or a point pattern with
  factor-valued marks, the levels of the factor are re-ordered so that
  the level \code{ref} is first and the others are moved down.
}
\usage{
  \method{relevel}{im}(x, ref, \dots)

  \method{relevel}{ppp}(x, ref, \dots)

  \method{relevel}{ppx}(x, ref, \dots)
}
\arguments{
  \item{x}{
    A pixel image (object of class \code{"im"}) with factor values,
    or a point pattern (object of class \code{"ppp"}, \code{"ppx"},
    \code{"lpp"} or \code{"pp3"}) with
    factor-valued marks.
  }
  \item{ref}{
    The reference level.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  These functions are methods for the generic
  \code{\link[stats]{relevel}}.

  If \code{x} is a pixel image (object of class \code{"im"}) with
  factor values, or a point pattern (object of class \code{"ppp"},
  \code{"ppx"}, \code{"lpp"} or \code{"pp3"})
  with factor-valued marks, the levels of the factor are changed so
  that the level specified by \code{ref} comes first. 
}
\value{
  Object of the same kind as \code{x}.
}
\author{
Adrian Baddeley
\email{Adrian.Baddeley@uwa.edu.au}
\url{http://www.maths.uwa.edu.au/~adrian/}
Rolf Turner
\email{r.turner@auckland.ac.nz}
and Ege Rubak
\email{rubak@math.aau.dk}
}
\seealso{
  \code{\link{mergeLevels}}
}
\examples{
  amacrine
  relevel(amacrine, "on")
}
\keyword{manip}
\keyword{spatial}
