\name{lixellate}
\alias{lixellate}
\title{
  Subdivide Segments of a Network
}
\description{
  Each line segment of a linear network
  will be divided into several shorter segments
  (line elements or lixels).
}
\usage{
lixellate(X, \dots, nsplit, eps, sparse = TRUE)
}
\arguments{
  \item{X}{
    A linear network (object of class \code{"linnet"})
    or a point pattern on a linear network (object of class \code{"lpp"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{nsplit}{
    Number of pieces into which \emph{each} line segment of \code{X}
    should be divided. Either a single integer, or an integer vector
    with one entry for each line segment in \code{X}.
    Incompatible with \code{eps}.
  }
  \item{eps}{
    Maximum length of the resulting pieces of line segment.
    A single numeric value.
    Incompatible with \code{nsplit}.
  }
  \item{sparse}{
    Optional. Logical value specifying whether the resulting
    linear network should be represented using a sparse matrix.
    If \code{sparse=NULL}, then the representation will be the
    same as in \code{X}.
  }
}
\details{
  Each line segment in \code{X} will be subdivided into equal pieces.
  The result is an object of the same kind as \code{X}, representing the
  same data as \code{X} except that the segments have been subdivided.

  Splitting is controlled by the arguments \code{nsplit} and \code{eps},
  exactly one of which should be given.

  If \code{nsplit} is given, it specifies the 
  number of pieces into which \emph{each} line segment of \code{X}
  should be divided. It should be either a single integer, or an integer vector
  of length equal to the number of line segments in \code{X}.

  If \code{eps} is given, it specifies the maximum length of
  any resulting piece of line segment. 

  It is strongly advisable to use \code{sparse=TRUE} (the default)
  to limit the computation time.

  If \code{X} is a point pattern (class \code{"lpp"}) then the
  spatial coordinates and marks of each data point are unchanged, but the
  local coordinates will change, because they are
  adjusted to map them to the new subdivided network.
}
\value{
  Object of the same kind as \code{X}.
}
\author{
  Greg McSwiggan, 
  \spatstatAuthors.
}
\seealso{
  \code{\link{linnet}}, \code{\link{lpp}}.
}
\examples{
   A <- lixellate(simplenet, nsplit=4)
   plot(A, main="lixellate(simplenet, nsplit=4)")
   points(vertices(A), pch=16)

   spiders
   lixellate(spiders, nsplit=3)
}
\keyword{spatial}
\keyword{manip}
