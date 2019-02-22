\name{venn.tess}
\alias{venn.tess}
\title{
  Tessellation Delimited by Several Sets
}
\description{
  Given a list of windows, construct the tessellation
  formed by all combinations of inclusion/exclusion of these
  windows. 
}
\usage{
venn.tess(\dots, window = NULL)
}
\arguments{
  \item{\dots}{
    Sets which delimit the tessellation.
    Any number of windows (objects of class \code{"owin"})
    or tessellations (objects of class \code{"tess"}).
  }
  \item{window}{
    Optional. The bounding window of the resulting tessellation.
    If not specified, the default is the union of all the
    arguments \code{\dots}.
  }
}
\details{
  The arguments \code{\dots} may be any number of windows.
  This function constructs a tessellation, like a Venn diagram,
  whose boundaries are made up of the boundaries of these sets.
  Each tile of the tessellation is defined by one of the 
  possible combinations in which each set is either included or
  excluded. 

  If the arguments \code{\dots} are named, then the resulting tiles
  will also have names, which identify the inclusion/exclusion
  combinations defining each tile. See the Examples.
}
\value{
  A tessellation (object of class \code{"tess"}).
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{intersect.tess}}.

  To construct other kinds of tessellations, see
  \code{\link{tess}},
  \code{\link{quadrats}}, \code{\link{hextess}},
  \code{\link{polartess}},
  \code{\link{dirichlet}}, \code{\link{delaunay}}, \code{\link{quantess}}
  and \code{\link{rpoislinetess}}.
}
\examples{
   V <- venn.tess(A=square(1),
                  B=square(c(-0.5, 0.5)),
                  window=square(c(-1,1.5)))
   V
   plot(V, do.labels=TRUE)
}
\keyword{spatial}
\keyword{math}
