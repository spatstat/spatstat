\name{gpc2owin}
\alias{gpc2owin}
\alias{owin2gpc}
\title{Convert Polygonal Region into Different Format} 
\description{
  Conversion between the different 
  representations of a polygonal region
  in the packages \pkg{spatstat} and \pkg{gpclib}.
}
\usage{
 gpc2owin(x)
 owin2gpc(x)
}
\arguments{
  \item{x}{Object representing a polygonal region.
    An object of class \code{"owin"} in the \pkg{spatstat}
    package (for \code{owin2gpc}) or an object of class
    \code{"gpc.poly"} in the \pkg{gpclib} package (for \code{gpc2owin}).
  }
}
\details{
  The packages \pkg{spatstat} and \pkg{gpclib} have slightly different
  internal formats for representing a polygonal region in the
  two-dimensional plane. In \pkg{gpclib} a polygonal region is
  an object of class \code{"gpc.poly"}, while in \pkg{spatstat}
  it is an object of class \code{"owin"} and of type \code{"polygonal"}.
  
  These two functions convert the two formats: \code{owin2gpc} converts
  an \code{"owin"} to a \code{"gpc.poly"}, while \code{gpc2owin}
  does the reverse.
  
  Conversion of a \code{"gpc.poly"} to an \code{"owin"} can also be
  performed by calling \code{\link{as.owin}}.
}
\value{
    An object of class \code{"owin"} in the \pkg{spatstat}
    package (for \code{gpc2owin}) or an object of class
    \code{"gpc.poly"} in the \pkg{gpclib} package (for \code{owin2gpc}).
}
\seealso{
  \code{\link{as.owin}}
}
\examples{
   data(letterR)
   if(spatstat.options("gpclib") && require(gpclib)) {
     R <- owin2gpc(letterR)
     L <- gpc2owin(R)
   } else cat("gpclib is not available\n")
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{utilities}
