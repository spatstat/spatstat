\name{Poisson}
\alias{Poisson}
\title{Poisson Point Process Model}
\description{
Creates an instance of the Poisson point process model
which can then be fitted to point pattern data.
}
\usage{
 Poisson()
}
\value{
  An object of class \code{"interact"}
  describing the interpoint interaction 
  structure of the Poisson point process
  (namely, there are no interactions).
}
\details{
  The function \code{\link{ppm}}, which fits point process models to 
  point pattern data, requires an argument \code{interaction}
  of class \code{"interact"}
  describing the interpoint interaction structure
  of the model to be fitted. 
  The appropriate description of the Poisson process is 
  provided by the value of the function \code{Poisson}.

  This works for all types of Poisson processes including
  multitype and nonstationary Poisson processes.
}
\seealso{
  \code{\link{ppm}},
  \code{\link{Strauss}}
}
\examples{
 data(nztrees)
 ppm(nztrees, ~1, Poisson())
 # fit the stationary Poisson process to 'nztrees'
 # no edge correction needed

 data(longleaf)
 \testonly{
   longleaf <- longleaf[seq(1, longleaf$n, by=50)]
 }
 longadult <- longleaf[longleaf$marks >= 30, ]
 longadult <- unmark(longadult)
 ppm(longadult, ~ x, Poisson())
 # fit the nonstationary Poisson process 
 # with intensity lambda(x,y) = exp( a + bx)

 data(lansing)
 # trees marked by species
 \testonly{
     lansing <- lansing[seq(1,lansing$n, by=30), ]
 }
 ppm(lansing, ~ marks, Poisson())
 # fit stationary marked Poisson process
 # with different intensity for each species

\dontrun{
 ppm(lansing, ~ marks * polynom(x,y,3), Poisson())
}
 # fit nonstationary marked Poisson process
 # with different log-cubic trend for each species
\testonly{
   # equivalent functionality - smaller dataset
   ppm(amacrine, ~ marks * polynom(x,y,2), Poisson())
}
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
 
 
