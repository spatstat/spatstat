\name{solist}
\alias{solist}
\title{
  List of Two-Dimensional Spatial Objects
}
\description{
  Make a list of two-dimensional spatial objects.
}
\usage{
solist(\dots, check = TRUE, promote = TRUE)
}
\arguments{
  \item{\dots}{
    Any number of objects, each representing a two-dimensional
    spatial dataset.
  }
  \item{check}{
    Logical value. If \code{TRUE}, check that each of the
    objects is a 2D spatial object.
  }
  \item{promote}{
    Logical value. If \code{TRUE}, test whether all objects belong to
    the \emph{same} class, and if so, promote the list of objects
    to the appropriate class of list.
  }
}
\details{
  This command creates an object of class \code{"solist"}
  (spatial object list)
  which represents a list of two-dimensional spatial datasets.
  The datasets do not necessarily belong to the same class.

  Typically the intention is that the datasets in the list
  should be treated in the same way, for example, they should
  be plotted side-by-side. The \pkg{spatstat} package
  provides a plotting function, \code{\link{plot.solist}},
  and many other functions for this class.
  
  In the \pkg{spatstat} package, various functions produce
  an object of class \code{"solist"}. For example, when
  a point pattern is split into several point patterns by
  \code{\link{split.ppp}}, or an image is split into several
  images by \code{\link{split.im}}, the result is of
  class \code{"solist"}.

  If \code{check=TRUE} then the code will check whether all
  objects belong to one of the classes of two-dimensional spatial
  objects defined in the \pkg{spatstat} package. If this is not true
  then the result will belong to the more general class \code{"anylist"}
  instead of \code{"solist"}. Set \code{check=FALSE}
  for efficiency, but only if you are sure that all the objects are valid.

  If \code{promote=TRUE} then the code will check whether all
  the objects \code{\dots} belong to the same class.
  If they are all point patterns (class \code{"ppp"}),
  the result will also belong to the class \code{"ppplist"}.
  If they are all pixel images (class \code{"im"}), the result
  will also belong to the class \code{"imlist"}.

  Use \code{\link{as.solist}} to convert a list to a \code{"solist"}.
}
\value{
  A list, usually belonging to the class \code{"solist"}.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\seealso{
  \code{\link{as.solist}}
}
\examples{
  solist(cells, density(cells))
  solist(cells, japanesepines, redwood)
}
\keyword{spatial}
\keyword{list}
\keyword{manip}
