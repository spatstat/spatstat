\name{is.rectangle}
\alias{is.rectangle}
\alias{is.polygonal}
\alias{is.mask}
\title{Determine Type of Window}
\description{
  Determine whether a window is a
  rectangle, a polygonal region, or a binary mask.
}
\usage{
  is.rectangle(w)
  is.polygonal(w)
  is.mask(w)
}
\arguments{
  \item{w}{
    Window to be inspected. An object of class \code{"owin"}.
  }
}
\value{
  Logical value, equal to \code{TRUE} if \code{w} is a window
  of the specified type.
}
\details{
  These simple functions determine whether a window \code{w}
  (object of class \code{"owin"}) is
  a rectangle (\code{is.rectangle(w) = TRUE}), 
  a domain with polygonal boundary (\code{is.polygonal(w) = TRUE}),
  or 
  a binary pixel mask (\code{is.mask(w) = TRUE}). 
}
\seealso{
  \code{\link{owin}}
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
