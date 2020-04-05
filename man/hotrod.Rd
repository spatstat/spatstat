\name{hotrod}
\alias{hotrod}
\title{
  Heat Kernel for a One-Dimensional Rod
}
\description{
  Calculate values of the heat kernel on a one-dimensional rod.
  The ends of the rod may be assumed to be insulated,
  or absorbing.
}
\usage{
hotrod(len, xsource, xquery, sigma, ends=c("insulated", "absorbing"), nmax=20)
}
\arguments{
  \item{len}{
    Length of the rod. A single number or numeric vector.
  }
  \item{xsource}{
    Positions of the source points, from the left end of the rod
    (in the same distance units as \code{len}).
    A single number or numeric vector.
  }
  \item{xquery}{
    Positions of the query points, from the left end of the rod
    (in the same distance units as \code{len}).
    A single number or numeric vector.
  }
  \item{sigma}{
    Bandwidth for kernel. A single number or a numeric vector.
  }
  \item{ends}{
    Character string (partially matched) specifying whether the
    ends of the rod are assumed to be insulated or absorbing.
  }
  \item{nmax}{
    Number of terms in the infinite sum to use.
    A single integer or an integer vector.
  }
}
\details{
  Computes the heat kernel as an infinite sum.
}
\value{
  Number or numeric vector.
}
\author{
  Greg McSwiggan and \adrian.
}
\examples{
  curve(hotrod(1, 0.1, x, 0.7))

  # check it's a probability density
  f <- function(x) hotrod(1, 0.1, x, 0.7)
  integrate(f, 0, 1)

  ## absorbing ends
  curve(hotrod(1, 0.1, x, 0.7, ends="a"))
}
\keyword{math}
