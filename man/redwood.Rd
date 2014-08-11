\name{redwood}
\alias{redwood}
\docType{data}
\title{
  California Redwoods Point Pattern (Ripley's Subset)
}
\description{
  Locations of 62 seedlings and saplings 
  of California redwood trees.
 
  The data represent the locations of 62 seedlings and saplings 
  of California redwood trees in a square sampling region.  
  They originate from Strauss (1975);
  the present data are a subset extracted by Ripley (1977) 
  in a subregion that has been rescaled to a unit square.

  The coordinates are rounded to the nearest 0.01 units,
  except for one point which has an \eqn{x} coordinate
  of 0.999, presumably to ensure that it is properly inside the window.
 
  There are many further analyses of this dataset. It is often used as a 
  canonical example of a clustered point pattern
  (see e.g. Diggle, 1983).

  The original, full redwood dataset is supplied in the \code{spatstat}
  library as \code{redwoodfull}.
} 
\format{
  An object of class \code{"ppp"}
  representing the point pattern of tree locations.
  The window has been rescaled to the unit square.

  See \code{\link{ppp.object}} for details of the format of a
  point pattern object.
}
\usage{data(redwood)}
\source{Strauss (1975), subset extracted by Ripley (1977)}
\seealso{
  \code{\link{redwoodfull}}
}
\references{
  Diggle, P.J. (1983)
  \emph{Statistical analysis of spatial point patterns}.
  Academic Press.

  Ripley, B.D. (1977)
  Modelling spatial patterns (with discussion).
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{39}, 172--212.

  Strauss, D.J. (1975)
  A model for clustering.
  \emph{Biometrika} \bold{63}, 467--475.
}
\keyword{datasets}
\keyword{spatial}
