\name{MultiHard}
\alias{MultiHard}
\title{The Multitype Hard Core Point Process Model}
\description{
Creates an instance of the multitype hard core point process model
which can then be fitted to point pattern data.
}
\usage{
  MultiHard(hradii, types=NULL)
}
\arguments{
  \item{hradii}{Matrix of hard core radii}
  \item{types}{Optional; vector of all possible types (i.e. the possible levels
    of the \code{marks} variable in the data)}
}
\value{
  An object of class \code{"interact"}
  describing the interpoint interaction
  structure of the multitype hard core process with
  hard core radii \eqn{hradii[i,j]}.
}
\details{
  This is a multitype version of the hard core process.
  A pair of points
  of types \eqn{i} and \eqn{j}
  must not lie closer than \eqn{h_{ij}}{h[i,j]} units apart.

  The argument \code{types} need not be specified in normal use.
  It will be determined automatically from the point pattern data set
  to which the MultiStrauss interaction is applied,
  when the user calls \code{\link{ppm}}. 
  However, the user should be confident that
  the ordering of types in the dataset corresponds to the ordering of
  rows and columns in the matrix \code{hradii}.

  The matrix \code{hradii}
  must be symmetric, with entries
  which are either positive numbers or \code{NA}. 
  A value of \code{NA} indicates that no distance constraint should be applied
  for this combination of types.
  
  Note that only the hardcore radii
  are specified in \code{MultiHard}.
  The canonical parameters \eqn{\log(\beta_j)}{log(beta[j])}
  are estimated by \code{\link{ppm}()}, not fixed in
  \code{MultiHard()}. 
}
\seealso{
  \code{\link{ppm}},
  \code{\link{pairwise.family}},
  \code{\link{ppm.object}},
  \code{\link{MultiStrauss}},
  \code{\link{MultiStraussHard}},
  \code{\link{Strauss}}
}
\examples{
   h <- matrix(c(1,2,2,1), nrow=2,ncol=2)

   # prints a sensible description of itself
   MultiHard(h)

   # Fit the stationary multitype hardcore process to `amacrine'
   # with hard core operating only between cells of the same type.
   h <- 0.02 * matrix(c(1, NA, NA, 1), nrow=2,ncol=2)
   ppm(amacrine ~1, MultiHard(h))
}
\section{Warnings}{
  In order that \code{\link{ppm}} can fit the multitype hard core
  model correctly to a point pattern \code{X}, this pattern must
  be marked, with \code{markformat} equal to \code{vector} and the
  mark vector \code{marks(X)} must be a factor.  If the argument
  \code{types} is specified it is interpreted as a set of factor
  levels and this set must equal \code{levels(marks(X))}.
}
\section{Changed Syntax}{
  Before \pkg{spatstat} version \code{1.37-0},
  the syntax of this function was different:
  \code{MultiHard(types=NULL, hradii)}.
  The new code attempts to handle the old syntax as well.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/},
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
}
\keyword{spatial}
\keyword{models}
