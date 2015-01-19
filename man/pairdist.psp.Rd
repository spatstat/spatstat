\name{pairdist.psp}
\alias{pairdist.psp}
\title{Pairwise distances between line segments}
\description{
  Computes the matrix of distances between all pairs of line segments
  in a line segment pattern.
}
\usage{
  \method{pairdist}{psp}(X, \dots, method="C", type="Hausdorff")
}
\arguments{
  \item{X}{
    A line segment pattern (object of class \code{"psp"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{method}{
    String specifying which method of calculation to use.
    Values are \code{"C"} and \code{"interpreted"}.
    Usually not specified.
  }
  \item{type}{
    Type of distance to be computed. Options are
    \code{"Hausdorff"} and \code{"separation"}. Partial matching is used.
  }
}
\value{
  A square matrix whose \code{[i,j]} entry is the distance
  between the line segments numbered \code{i} and \code{j}.
}
\details{
  This function computes the distance between each pair of
  line segments in \code{X}, and returns the matrix of distances.
  
  This is a method for the generic function \code{\link{pairdist}}
  for the class \code{"psp"}.

  The distances between line segments are measured in one of two ways:
  \itemize{
    \item if \code{type="Hausdorff"}, distances are computed
    in the Hausdorff metric. The Hausdorff
    distance between two line segments is the \emph{maximum} distance
    from any point on one of the segments to the nearest point on
    the other segment.
    \item if \code{type="separation"}, distances are computed
    as the \emph{minimum} distance from a point on one line segment to
    a point on the other line segment. For example, line segments which
    cross over each other have separation zero.
  }
  
  The argument \code{method} is not normally used. It is
  retained only for checking the validity of the software.
  If \code{method = "interpreted"} then the distances are
  computed using interpreted \R code only. If \code{method="C"}
  (the default) then compiled \code{C} code is used,
  which is somewhat faster.
}
\seealso{
  \code{\link{crossdist}},
  \code{\link{nndist}},
  \code{\link{pairdist.ppp}}
}
\examples{
   L <- psp(runif(10), runif(10), runif(10), runif(10), owin())
   D <- pairdist(L)
   S <- pairdist(L, type="sep")
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}

