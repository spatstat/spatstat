\name{unstack.msr}
\alias{unstack.msr}
\title{
  Separate a Vector Measure into its Scalar Components
}
\description{
  Converts a vector-valued measure into a list of
  scalar-valued measures. 
}
\usage{
\method{unstack}{msr}(x, \dots)
}
\arguments{
  \item{x}{
    A measure (object of class \code{"msr"}).
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  This is a method for the generic
  \code{\link[utils]{unstack}} for the class \code{"msr"} of measures. 
  
  If \code{x} is a vector-valued measure, then 
  \code{y <- unstack(x)} is a list of scalar-valued measures
  defined by the components of \code{x}. 
  The \code{j}th entry of the list, \code{y[[j]]}, is equivalent to
  the \code{j}th component of the vector measure \code{x}.
  
  If \code{x} is a scalar-valued measure, then
  the result is a list consisting of one entry, which is \code{x}.
}
\value{
  A list of measures, of class \code{"solist"}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link[utils]{unstack}}

  \code{\link{unstack.ppp}}

  \code{\link{split.msr}}.
}
\examples{
   fit <- ppm(cells ~ x)
   m <- residuals(fit, type="score")
   m
   unstack(m)
}
\keyword{spatial}
\keyword{manip}
