\name{quad.object}
\alias{quad.object} %DoNotExport
\title{Class of Quadrature Schemes}
\description{
  A class \code{"quad"} to represent a quadrature scheme.
}
\details{
  A (finite) quadrature scheme is a list of quadrature points
  \eqn{u_j}{u[j]} and associated weights \eqn{w_j}{w[j]}
  which is used to approximate an integral by a finite sum:
  \deqn{
    \int f(x) dx \approx \sum_j f(u_j) w_j
  }{
    integral(f(x) dx) ~= sum( f(u[j]) w[j] )
  }
  Given a point pattern dataset, a \emph{Berman-Turner}
  quadrature scheme is one which includes all these data points,
  as well as a nonzero number of other (``dummy'') points.

  These quadrature schemes are used to approximate the
  pseudolikelihood of a point process, in the method of
  Baddeley and Turner (2000) (see Berman and Turner (1992)).
  Accuracy and computation time both increase with the number of
  points in the quadrature scheme.

  An object of class \code{"quad"} represents a Berman-Turner
  quadrature scheme. It can be passed as an argument to
  the model-fitting function \code{\link{ppm}}, which
  requires a quadrature scheme.

  An object of this class contains at least the following elements:
  \tabular{ll}{
    \code{data}: \tab an object of class \code{"ppp"} \cr
                 \tab giving the locations (and marks) of the data points.\cr
    \code{dummy}: \tab an object of class \code{"ppp"} \cr
                 \tab giving the locations (and marks) of the dummy points.\cr
    \code{w}: \tab vector of nonnegative weights for the quadrature points\cr
  }
  Users are strongly advised not to manipulate these entries
  directly.
  
  The domain of quadrature is specified by \code{Window(dummy)}
  while the observation window (if this needs to be specified
  separately) is taken to be \code{Window(data)}. 

  The weights vector \code{w} may also have an attribute
  \code{attr(w, "zeroes")} equivalent to the logical vector
  \code{(w == 0)}. If this is absent then all points are known to
  have positive weights.
  
  To create an object of class \code{"quad"},
  users would typically call the
  high level function \code{\link{quadscheme}}.
  (They are actually
  created by the low level function \code{quad}.)
  
  Entries are extracted from a \code{"quad"} object by the functions
  \code{x.quad},
  \code{y.quad},
  \code{w.quad} and
  \code{marks.quad},
  which extract the \eqn{x} coordinates, \eqn{y} coordinates,
  weights, and marks, respectively. The function
  \code{n.quad} returns the total number of quadrature points
  (dummy plus data).

  An object of class \code{"quad"} can be converted into an ordinary
  point pattern by the function \code{\link{union.quad}} which simply
  takes the union of the data and dummy points.

  Quadrature schemes can be plotted using \code{\link{plot.quad}}
  (a method for the generic \code{\link{plot}}).
}
\seealso{
  \code{\link{quadscheme}},
  \code{\link{ppm}}
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{attribute}
 
 
