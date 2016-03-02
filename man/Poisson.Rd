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
 ppm(nztrees ~1, Poisson())
 # fit the stationary Poisson process to 'nztrees'
 # no edge correction needed

 lon <- longleaf
 \testonly{
   lon <- lon[seq(1, npoints(lon), by=50)]
 }
 longadult <- unmark(subset(lon, marks >= 30))
 ppm(longadult ~ x, Poisson())
 # fit the nonstationary Poisson process 
 # with intensity lambda(x,y) = exp( a + bx)

 # trees marked by species
 lans <- lansing
 \testonly{
     lans <- lans[seq(1, npoints(lans), by=30)]
 }
 ppm(lans ~ marks, Poisson())
 # fit stationary marked Poisson process
 # with different intensity for each species

\dontrun{
 ppm(lansing ~ marks * polynom(x,y,3), Poisson())
}
 # fit nonstationary marked Poisson process
 # with different log-cubic trend for each species
\testonly{
   # equivalent functionality - smaller dataset
   ppm(amacrine ~ marks * polynom(x,y,2), Poisson())
}
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au},
 Rolf Turner \email{r.turner@auckland.ac.nz}
 and Ege Rubak \email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{models}
 
 
