\name{delaunayNetwork}
\alias{delaunayNetwork}
\alias{dirichletNetwork}
\title{
  Linear Network of Delaunay Triangulation or Dirichlet Tessellation
}
\description{
  Computes the edges of the Delaunay triangulation
  or Dirichlet tessellation of a point pattern,
  and returns the result as a linear network object.
}
\usage{
delaunayNetwork(X)

dirichletNetwork(X, \dots)
}
\arguments{
  \item{X}{A point pattern (object of class \code{"ppp"}).}
  \item{\dots}{Arguments passed to \code{\link{as.linnet.psp}}}
}
\details{
  For \code{delaunayNetwork}, points of \code{X} which are neighbours
  in the Delaunay triangulation
  (see \code{\link{delaunay}}) will be joined by a straight line.
  The result will be returned as a linear network (object of class
  \code{"linnet"}).

  For \code{dirichletNetwork}, the Dirichlet tessellation is computed
  (see \code{\link{dirichlet}}) and the edges of the
  tiles of the tessellation are extracted. This is converted to a linear
  network using \code{\link{as.linnet.psp}}.
}
\value{
  Linear network (object of class \code{"linnet"})
  or \code{NULL}.
}
\author{
Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}


Rolf Turner \email{r.turner@auckland.ac.nz}

and Ege Rubak \email{rubak@math.aau.dk}

}
\seealso{
  \code{\link{delaunay}},
  \code{\link{dirichlet}},
  \code{\link{delaunayDistance}}
}
\examples{
  LE <- delaunayNetwork(cells)
  LI <- dirichletNetwork(cells)
}
\keyword{spatial}
\keyword{manip}
