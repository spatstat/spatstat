\name{bw.stoyan}
\alias{bw.stoyan}
\title{
  Stoyan's Rule of Thumb for Bandwidth Selection
}
\description{
  Computes a rough estimate of the appropriate bandwidth
  for kernel smoothing estimators of the pair correlation function
  and other quantities.
}
\usage{
bw.stoyan(X, co=0.15)
}

\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
  \item{co}{
    Coefficient appearing in the rule of thumb. See Details.
  }
}
\details{
  Estimation of the pair correlation function and other quantities
  by smoothing methods requires a choice of the smoothing bandwidth.
  Stoyan and Stoyan (1995, equation (15.16), page 285) proposed a
  rule of thumb for choosing the smoothing bandwidth.
  
  For the Epanechnikov kernel, the rule of thumb is to set
  the kernel's half-width \eqn{h} to
  \eqn{0.15/\sqrt{\lambda}}{0.15/sqrt(lambda)} where
  \eqn{\lambda}{lambda} is the estimated intensity of the point pattern,
  typically computed as the number of points of \code{X} divided by the
  area of the window containing \code{X}.

  For a general kernel, the corresponding rule is to set the
  standard deviation of the kernel to
  \eqn{\sigma = 0.15/\sqrt{5\lambda}}{sigma = 0.15/sqrt(5 * lambda)}.
  
  The coefficient \eqn{0.15} can be tweaked using the
  argument \code{co}. 
}
\value{
  A numerical value giving the selected bandwidth (the standard
  deviation of the smoothing kernel).
}
\seealso{
  \code{\link{pcf}},
  \code{\link{bw.relrisk}}
}
\examples{
  data(shapley)
  bw.stoyan(shapley)
}
\references{
  Stoyan, D. and Stoyan, H. (1995)
  Fractals, random shapes and point fields:
  methods of geometrical statistics.
  John Wiley and Sons.
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
