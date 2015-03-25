\name{marktable}
\alias{marktable}
\title{Tabulate Marks in Neighbourhood of Every Point in a Point Pattern}
\description{
  Visit each point in a point pattern, find the neighbouring points,
  and compile a frequency table of the marks of these neighbour points.
}
\usage{
   marktable(X, R, N, exclude=TRUE, collapse=FALSE)
}
\arguments{
  \item{X}{
    A marked point pattern. 
    An object of class \code{"ppp"}.
  }
  \item{R}{
    Neighbourhood radius. Incompatible with \code{N}.
  }
  \item{N}{
    Number of neighbours of each point. Incompatible with \code{R}.
  }
  \item{exclude}{
    Logical. If \code{exclude=TRUE}, the neighbours of a point
    do not include the point itself. If \code{exclude=FALSE},
    a point belongs to its own neighbourhood.
  }
  \item{collapse}{
    Logical. If \code{collapse=FALSE} (the default) the results for
    each point are returned as separate rows of a table.
    If \code{collapse=TRUE}, the results are aggregated according to the
    type of point.
  }
}
\value{
  A contingency table (object of class \code{"table"}).
  If \code{collapse=FALSE}, the table has one row for
  each point in \code{X}, and one column for each possible mark value.
  If \code{collapse=TRUE}, the table has one row and one column
  for each possible mark value.
}
\details{
  This algorithm visits each point in the point pattern \code{X},
  inspects all the neighbouring points within a radius \code{R} of the current
  point (or the \code{N} nearest neighbours of the current point),
  and compiles a frequency table of the marks attached to the
  neighbours. 

  The dataset \code{X} must be a multitype point pattern, that is,
  \code{marks(X)} must be a \code{factor}.
  
  If \code{collapse=FALSE} (the default), 
  the result is a two-dimensional contingency table with one row for
  each point in the pattern, and one column for each possible mark
  value. The \code{[i,j]} entry in the table gives the number of
  neighbours of point \code{i} that have mark \code{j}.

  If \code{collapse=TRUE}, this contingency table is aggregated
  according to the type of point, so that the result is a contingency
  table with one row and one column for each possible mark value.
  The  \code{[i,j]} entry in the table gives the number of
  neighbours of a point with mark \code{i} that have mark \code{j}.

  To perform more complicated calculations on the neighbours of every
  point, use \code{\link{markstat}} or \code{\link{applynbd}}. 
}
\seealso{
  \code{\link{markstat}},
  \code{\link{applynbd}},
  \code{\link{Kcross}},
  \code{\link{ppp.object}},
  \code{\link{table}}
}
\examples{
  head(marktable(amacrine, 0.1))
  head(marktable(amacrine, 0.1, exclude=FALSE))
  marktable(amacrine, N=1, collapse=TRUE)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{programming}
 
