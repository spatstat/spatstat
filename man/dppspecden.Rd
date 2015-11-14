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
  \adrian
  
  
  \rolf
  
  and \ege
  
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
