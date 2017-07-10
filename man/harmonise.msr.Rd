\name{harmonise.msr}
\alias{harmonise.msr}
\title{Make Measures Compatible}
\description{
  Convert several measures to a common quadrature scheme
}
\usage{
\method{harmonise}{msr}(\dots)
}
\arguments{
  \item{\dots}{
    Any number of measures (objects of class \code{"msr"}).
  } 
}
\details{
  This function makes any number of measures compatible,
  by converting them all to a common quadrature scheme.

  The command \code{\link{harmonise}} is generic. This is the
  method for objects of class \code{"msr"}.
}
\value{
  A list, of length equal to the number of arguments \code{\dots},
  whose entries are measures.
}
\author{
  \spatstatAuthors.
}
\examples{
  fit1 <- ppm(cells ~ x)
  fit2 <- ppm(rpoispp(ex=cells) ~ x)
  m1 <- residuals(fit1)
  m2 <- residuals(fit2)
  harmonise(m1, m2)
  s1 <- residuals(fit1, type="score")
  s2 <- residuals(fit2, type="score")
  harmonise(s1, s2)
}
\seealso{
  \code{\link{harmonise}},
  \code{\link{msr}}
}
\keyword{spatial}
\keyword{manip}
