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
  \adrian
  
  
  \rolf
  
  and \ege
  
}





