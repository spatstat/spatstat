\name{duplicated.ppp}
\alias{duplicated.ppp}
\alias{duplicated.ppx}
\alias{anyDuplicated.ppp}
\alias{anyDuplicated.ppx}
\title{Determine Duplicated Points in a Spatial Point Pattern}
\description{
  Determines which points in a spatial point pattern
  are duplicates of previous points, and returns a logical vector.
}
\usage{
 \method{duplicated}{ppp}(x, \dots, rule=c("spatstat", "deldir", "unmark"))

 \method{duplicated}{ppx}(x, \dots)

 \method{anyDuplicated}{ppp}(x, \dots)

 \method{anyDuplicated}{ppx}(x, \dots)
}
\arguments{
  \item{x}{
    A spatial point pattern
    (object of class \code{"ppp"} or \code{"ppx"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{rule}{
    Character string.
    The rule for determining duplicated points.
  }
}
\value{
  \code{duplicated(x)} returns
  a logical vector of length equal to the number of points in \code{x}.

  \code{anyDuplicated(x)} is a number equal to 0 if there are no
  duplicated points, and otherwise is equal to the index of the first
  duplicated point.
}
\details{
  These are methods for the generic functions \code{\link{duplicated}}
  and \code{\link{anyDuplicated}} for 
  point pattern datasets (of class \code{"ppp"}, see
  \code{\link{ppp.object}}, or class \code{"ppx"}).

  \code{anyDuplicated(x)} is a faster version of
  \code{any(duplicated(x))}.

  Two points in a point pattern are deemed to be identical
  if their \eqn{x,y} coordinates are the same,
  and their marks are also the same (if they carry marks).
  The Examples section illustrates how it is possible for
  a point pattern to contain a pair of identical points.

  This function determines which points in \code{x} duplicate
  other points that appeared earlier in the sequence. It 
  returns a logical vector with entries that are \code{TRUE}
  for duplicated points and \code{FALSE} for unique (non-duplicated)
  points.

  If \code{rule="spatstat"} (the default), two points are deemed
  identical if their coordinates are equal according to \code{==},
  \emph{and} their marks are equal according to \code{==}.
  This is the most stringent possible test.
  If \code{rule="unmark"}, duplicated points are
  determined by testing equality of their coordinates only,
  using \code{==}.
  If \code{rule="deldir"}, duplicated points are
  determined by testing equality of their coordinates only,
  using the function \code{\link[deldir]{duplicatedxy}}
  in the package \pkg{deldir}, which currently uses
  \code{\link{duplicated.data.frame}}. Setting \code{rule="deldir"}
  will ensure consistency with functions in the \pkg{deldir} package.
}
\seealso{
  \code{\link{ppp.object}},
  \code{\link{unique.ppp}},
  \code{\link{multiplicity.ppp}}
}
\examples{
   X <- ppp(c(1,1,0.5), c(2,2,1), window=square(3))
   duplicated(X)
   duplicated(X, rule="deldir")
}

\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{methods}
 
