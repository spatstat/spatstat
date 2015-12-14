\name{unique.ppp}
\alias{unique.ppp}
\alias{unique.ppx}
\title{Extract Unique Points from a Spatial Point Pattern}
\description{
  Removes any points that are identical to other points
  in a spatial point pattern.
}
\usage{
 \method{unique}{ppp}(x, \dots, warn=FALSE)

 \method{unique}{ppx}(x, \dots, warn=FALSE)
}
\arguments{
  \item{x}{
    A spatial point pattern
    (object of class \code{"ppp"} or \code{"ppx"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link{duplicated.ppp}}
    or \code{\link{duplicated.data.frame}}.
  }
  \item{warn}{
    Logical. If \code{TRUE}, issue a warning message if any
    duplicated points were found.
  }
}
\value{
  Another point pattern object.
}
\details{
  These are methods for the generic function \code{unique} for 
  point pattern datasets (of class \code{"ppp"}, see
  \code{\link{ppp.object}}, or class \code{"ppx"}).

  This function removes duplicate points in \code{x},
  and returns a point pattern.

  Two points in a point pattern are deemed to be identical
  if their \eqn{x,y} coordinates are the same,
  \emph{and} their marks are the same (if they carry marks).
  This is the default rule: see \code{\link{duplicated.ppp}}
  for other options. 
}
\seealso{
  \code{\link{ppp.object}},
  \code{\link{duplicated.ppp}},
  \code{\link{multiplicity.ppp}}
}
\examples{
   X <- ppp(c(1,1,0.5), c(2,2,1), window=square(3))
   unique(X)
   unique(X, rule="deldir")
}

\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{methods}
 
