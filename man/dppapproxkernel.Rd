\name{dppapproxkernel}
\alias{dppapproxkernel}
\title{Approximate Determinantal Point Process Kernel}
\description{
  Returns an approximation to the kernel of a determinantal
  point process, as a function of one argument \eqn{x}.
}
\usage{dppapproxkernel(model, trunc = 0.99, W = NULL)}
\arguments{
  \item{model}{Object of class \code{"detpointprocfamily"}.}
  \item{trunc}{Numeric specifying how the model truncation is
    performed. See Details section of \code{\link{simulate.detpointprocfamily}}.
  }
  \item{W}{Optional window -- undocumented at the moment.}
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





