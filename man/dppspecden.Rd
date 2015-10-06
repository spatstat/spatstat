\name{dppspecden}
\alias{dppspecden}
\title{Extract Spectral Density from Determinantal Point Process Model Object}
\description{
  Returns the spectral density of a determinantal point process
  model as a function of one argument \code{x}.
}
\usage{dppspecden(model)}
\arguments{
  \item{model}{Model of class \code{"detpointprocfamily"}.}
}
\value{A function}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\examples{
model <- dppMatern(lambda = 100, alpha=.01, nu=1, d=2)
dppspecden(model)
}
\seealso{
  \code{\link{dppspecdenrange}}
}
\keyword{spatial}
\keyword{models}
