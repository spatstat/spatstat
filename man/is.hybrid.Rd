\name{is.hybrid}  %DontDeclareMethods
\alias{is.hybrid}
\alias{is.hybrid.ppm}
\alias{is.hybrid.interact}
\title{
  Test Whether Object is a Hybrid
}
\description{
  Tests where a point process model or point process interaction
  is a hybrid of several interactions.
}
\usage{
is.hybrid(x)

\method{is.hybrid}{ppm}(x)

\method{is.hybrid}{interact}(x)
}
\arguments{
  \item{x}{
    A point process model (object of class \code{"ppm"})
    or a point process interaction structure
    (object of class \code{"interact"}).
  }
}
\details{
  A \emph{hybrid} (Baddeley, Turner, Mateu and Bevan, 2012)
  is a point process model created by combining two or more
  point process models, or an interpoint interaction created by combining
  two or more interpoint interactions.
  
  The function \code{is.hybrid} is generic, with methods for
  point process models (objects of class \code{"ppm"})
  and point process interactions
  (objects of class \code{"interact"}).
  These functions return \code{TRUE} if the object \code{x} is a hybrid,
  and \code{FALSE} if it is not a hybrid.

  Hybrids of two or more interpoint interactions
  are created by the function \code{\link{Hybrid}}.
  Such a hybrid interaction can then be fitted to point pattern data
  using \code{\link{ppm}}. 
}
\value{
  \code{TRUE} if the object is a hybrid, and \code{FALSE} otherwise.
}
\references{
  Baddeley, A., Turner, R., Mateu, J. and Bevan, A. (2013)
  Hybrids of Gibbs point process models and their implementation.
  \emph{Journal of Statistical Software} \bold{55}:11, 1--43.
  \url{http://www.jstatsoft.org/v55/i11/}
}
\seealso{
  \code{\link{Hybrid}}
}
\examples{
  S <- Strauss(0.1)
  is.hybrid(S)
  H <- Hybrid(Strauss(0.1), Geyer(0.2, 3))
  is.hybrid(H)

  data(redwood)
  fit <- ppm(redwood, ~1, H)
  is.hybrid(fit)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{models}
