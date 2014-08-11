\name{MultiStrauss}
\alias{MultiStrauss}
\title{The Multitype Strauss Point Process Model}
\description{
Creates an instance of the multitype Strauss point process model
which can then be fitted to point pattern data.
}
\usage{
  MultiStrauss(types=NULL, radii)
}
\arguments{
  \item{types}{Optional; vector of all possible types (i.e. the possible levels
    of the \code{marks} variable in the data)}
  \item{radii}{Matrix of interaction radii}
}
\value{
  An object of class \code{"interact"}
  describing the interpoint interaction
  structure of the multitype Strauss process with
  interaction radii \eqn{radii[i,j]}.
}
\details{
  The (stationary) multitype
  Strauss process with \eqn{m} types, with interaction radii
  \eqn{r_{ij}}{r[i,j]} and 
  parameters \eqn{\beta_j}{beta[j]} and \eqn{\gamma_{ij}}{gamma[i,j]}
  is the pairwise interaction point process
  in which each point of type \eqn{j}
  contributes a factor \eqn{\beta_j}{beta[j]} to the 
  probability density of the point pattern, and a pair of points
  of types \eqn{i} and \eqn{j} closer than \eqn{r_{ij}}{r[i,j]}
  units apart contributes a factor
  \eqn{\gamma_{ij}}{gamma[i,j]} to the density.

  The nonstationary multitype Strauss process is similar except that 
  the contribution of each individual point \eqn{x_i}{x[i]}
  is a function \eqn{\beta(x_i)}{beta(x[i])}
  of location and type, rather than a constant beta. 
 
  The function \code{\link{ppm}()}, which fits point process models to 
  point pattern data, requires an argument 
  of class \code{"interact"} describing the interpoint interaction
  structure of the model to be fitted. 
  The appropriate description of the multitype
  Strauss process pairwise interaction is
  yielded by the function \code{MultiStrauss()}. See the examples below.

  The argument \code{types} need not be specified in normal use.
  It will be determined automatically from the point pattern data set
  to which the MultiStrauss interaction is applied,
  when the user calls \code{\link{ppm}}. 
  However, the user should be confident that
  the ordering of types in the dataset corresponds to the ordering of
  rows and columns in the matrix \code{radii}.

  The matrix \code{radii} must be symmetric, with entries
  which are either positive numbers or \code{NA}. 
  A value of \code{NA} indicates that no interaction term should be included
  for this combination of types.
  
  Note that only the interaction radii are
  specified in \code{MultiStrauss}.  The canonical
  parameters \eqn{\log(\beta_j)}{log(beta[j])} and
  \eqn{\log(\gamma_{ij})}{log(gamma[i,j])} are estimated by
  \code{\link{ppm}()}, not fixed in \code{MultiStrauss()}.
}
\seealso{
  \code{\link{ppm}},
  \code{\link{pairwise.family}},
  \code{\link{ppm.object}},
  \code{\link{Strauss}},
  \code{\link{MultiHard}}
}
\examples{
   r <- matrix(c(1,2,2,1), nrow=2,ncol=2)
   MultiStrauss(radii=r)
   # prints a sensible description of itself
   r <- 0.03 * matrix(c(1,2,2,1), nrow=2,ncol=2)
   X <- amacrine
   \testonly{
      X <- X[ owin(c(0, 0.8), c(0, 1)) ]
   }
   ppm(X, ~1, MultiStrauss(, r))
   # fit the stationary multitype Strauss process to `amacrine'
   # Note the comma; needed since "types" is not specified.

   \dontrun{
   ppm(X, ~polynom(x,y,3), MultiStrauss(c("off","on"), r))
   # fit a nonstationary multitype Strauss process with log-cubic trend
   }
}
\section{Warnings}{
  In order that \code{\link{ppm}} can fit the multitype Strauss
  model correctly to a point pattern \code{X}, this pattern must
  be marked, with \code{markformat} equal to \code{vector} and the
  mark vector \code{marks(X)} must be a factor.  If the argument
  \code{types} is specified it is interpreted as a set of factor
  levels and this set must equal \code{levels(marks(X))}.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
