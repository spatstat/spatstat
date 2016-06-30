\name{CDF}
\alias{CDF}
\alias{CDF.density}
\title{
  Cumulative Distribution Function From Kernel Density Estimate
}
\description{
  Given a kernel estimate of a probability density,
  compute the corresponding cumulative distribution function.
}
\usage{
CDF(f, \dots)

\method{CDF}{density}(f, \dots, warn = TRUE)
}
\arguments{
  \item{f}{
    Density estimate (object of class \code{"density"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{warn}{
    Logical value indicating whether to issue a warning if the
    density estimate \code{f} had to be renormalised because it
    was computed in a restricted interval.
  }
}
\details{
  \code{CDF} is generic, with a method for class \code{"density"}.
  
  This calculates the cumulative distribution function
  whose probability density has been estimated and stored in the object
  \code{f}. The object \code{f} must belong to the class \code{"density"},
  and would typically have been obtained from a call to the function
  \code{\link[stats]{density}}.
}
\value{
  A function, which can be applied to any numeric value or vector of
  values.
}
\author{
  \spatstatAuthors
}
\seealso{
  \code{\link[stats]{density}},
  \code{\link{quantile.density}}
}
\examples{
   b <- density(runif(10))
   f <- CDF(b)
   f(0.5)
   plot(f)
}
\keyword{nonparametric}
\keyword{univar}
