\name{LambertW}
\alias{LambertW}
\title{
  Lambert's W Function
}
\description{
  Computes Lambert's W-function.
}
\usage{
LambertW(x)
}
\arguments{
  \item{x}{
    Vector of nonnegative numbers.
  }
}
\details{
  Lambert's W-function is the inverse function of
  \eqn{f(y) = y e^y}{f(y) = y * exp(y)}.
  That is, \eqn{W} is the function such that
  \deqn{
    W(x) e^{W(x)} = x
  }{
    W(x) * exp(W(x)) = x
  }

  This command \code{LambertW} computes \eqn{W(x)} for each entry
  in the argument \code{x}. 
  If the library \pkg{gsl} has been installed, then the function
  \code{lambert_W0} in that library is invoked. Otherwise,
  values of the W-function are computed by root-finding, using the
  function \code{\link[stats]{uniroot}}.

  Computation using \pkg{gsl} is about 100 times faster.

  If any entries of \code{x} are infinite or \code{NA}, the corresponding
  results are \code{NA}.
}
\value{
  Numeric vector.
}
\references{
  Corless, R, Gonnet, G, Hare, D, Jeffrey, D and Knuth, D (1996),
  On the Lambert W function.
  \emph{Computational Mathematics}, \bold{5}, 325--359.
  
  Roy, R and Olver, F (2010),
  Lambert W function. In Olver, F, Lozier, D and Boisvert, R (eds.),
  \emph{{NIST} Handbook of Mathematical Functions},
  Cambridge University Press.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\examples{
   LambertW(exp(1))
}
\keyword{math}
