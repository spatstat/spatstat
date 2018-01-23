\name{test.crossing.psp}
\alias{test.crossing.psp}
\alias{test.selfcrossing.psp}
\title{
  Check Whether Segments Cross
}
\description{
  Determine whether there is a crossing (intersection) between
  each pair of line segments.
}
\usage{
test.crossing.psp(A, B)
test.selfcrossing.psp(A)
}
\arguments{
  \item{A,B}{
    Line segment patterns (objects of class \code{"psp"}).
  }
}
\details{
  These functions decide whether the given line segments intersect
  each other.

  If \code{A} and \code{B} are two spatial patterns of line segments,
  \code{test.crossing.psp(A, B)}
  returns a logical matrix in which the entry on row \code{i}, column
  \code{j} is equal to \code{TRUE} if segment \code{A[i]}
  has an intersection with segment \code{B[j]}.

  If \code{A} is a pattern of line segments, \code{test.selfcross.psp(A)}
  returns a symmetric logical matrix in which the entry on row \code{i}, column
  \code{j} is equal to \code{TRUE} if segment \code{A[i]}
  has an intersection with segment \code{A[j]}.
}
\value{
  A logical matrix.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{psp}}
}
\examples{
  B <- edges(letterR)
  A <- rpoisline(5, Frame(B))
  MA <- test.selfcrossing.psp(A)
  MAB <- test.crossing.psp(A, B)
}
\keyword{spatial}
\keyword{math}
