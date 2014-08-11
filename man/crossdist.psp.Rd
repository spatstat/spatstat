\name{crossdist.psp}
\alias{crossdist.psp}
\title{Pairwise distances between two different line segment patterns}
\description{
  Computes the distances between all pairs of line segments
  taken from two different line segment patterns.
}
\usage{
  \method{crossdist}{psp}(X, Y, \dots, method="Fortran", type="Hausdorff")
}
\arguments{
  \item{X,Y}{
    Line segment patterns (objects of class \code{"psp"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{method}{String specifying which method of calculation to use.
    Values are \code{"Fortran"}, \code{"C"} and \code{"interpreted"}.
    Usually not specified.
  }
  \item{type}{
    Type of distance to be computed. Options are
    \code{"Hausdorff"} and \code{"separation"}. Partial matching is used.
  }
}
\value{
  A matrix whose \code{[i,j]} entry is the distance
  from the \code{i}-th line segment in \code{X}
  to the \code{j}-th line segment in \code{Y}.
}
\details{
  This is a method for the generic function \code{\link{crossdist}}.

  Given two line segment patterns, 
  this function computes the distance from each line segment
  in the first pattern to each line segment in the second pattern,
  and returns a matrix containing these distances.

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
  computed using interpreted R code only. If \code{method="Fortran"}
  (the default) then Fortran code is used. 
  The Fortran code is several times faster.
}
\seealso{
  \code{\link{pairdist}},
  \code{\link{nndist}},
  \code{\link{Gest}}
}
\examples{
   L1 <- psp(runif(5), runif(5), runif(5), runif(5), owin())
   L2 <- psp(runif(10), runif(10), runif(10), runif(10), owin())
   D <- crossdist(L1, L2)
   #result is a 5 x 10 matrix
   S <- crossdist(L1, L2, type="sep")
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}

\keyword{math}
