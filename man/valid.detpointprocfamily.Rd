\name{valid.detpointprocfamily}
\alias{valid.detpointprocfamily}
\title{Check Validity of a Determinantal Point Process Model}
\description{
  Checks the validity of a determinantal point process model.
}
\usage{
   \method{valid}{detpointprocfamily}(object, \dots)
}
\arguments{
  \item{object}{Model of class \code{"detpointprocfamily"}.}
  \item{\dots}{Ignored.}
}
\value{Logical}
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
model1 <- dppMatern(lambda=100, alpha=.01, nu=1, d=2)
valid(model1)
model2 <- dppMatern(lambda=100, alpha=1, nu=1, d=2)
valid(model2)
}
\seealso{
  \code{\link{valid}}
}
\keyword{spatial}
\keyword{models}
