\name{MultiStraussHard}
\alias{MultiStraussHard}
\title{The Multitype/Hard Core Strauss Point Process Model}
\description{
Creates an instance of the multitype/hard core Strauss point process model
which can then be fitted to point pattern data.
}
\usage{
  MultiStraussHard(iradii, hradii, types=NULL)
}
\arguments{
  \item{iradii}{Matrix of interaction radii}
  \item{hradii}{Matrix of hard core radii}
  \item{types}{Optional; vector of all possible types (i.e. the possible levels
    of the \code{marks} variable in the data)}
}
\value{
  An object of class \code{"interact"}
  describing the interpoint interaction
  structure of the multitype/hard core Strauss process with
  interaction radii \eqn{iradii[i,j]}
  and hard core radii \eqn{hradii[i,j]}.
}
\details{
  This is a hybrid of the multitype Strauss process
  (see \code{\link{MultiStrauss}}) and the hard core process
  (case \eqn{\gamma=0}{gamma = 0} of the Strauss process).
  A pair of points
  of types \eqn{i} and \eqn{j}
  must not lie closer than \eqn{h_{ij}}{h[i,j]} units apart;
  if the pair lies more than \eqn{h_{ij}}{h[i,j]} and less than 
  \eqn{r_{ij}}{r[i,j]} units apart, it
  contributes a factor
  \eqn{\gamma_{ij}}{gamma[i,j]} to the probability density.

  The argument \code{types} need not be specified in normal use.
  It will be determined automatically from the point pattern data set
  to which the MultiStraussHard interaction is applied,
  when the user calls \code{\link{ppm}}. 
  However, the user should be confident that
  the ordering of types in the dataset corresponds to the ordering of
  rows and columns in the matrices \code{iradii} and \code{hradii}.

  The matrices \code{iradii} and \code{hradii}
  must be symmetric, with entries
  which are either positive numbers or \code{NA}. 
  A value of \code{NA} indicates that no interaction term should be included
  for this combination of types.
  
  Note that only the interaction radii and hardcore radii
  are specified in \code{MultiStraussHard}.
  The canonical parameters \eqn{\log(\beta_j)}{log(beta[j])}
  and \eqn{\log(\gamma_{ij})}{log(gamma[i,j])}
  are estimated by \code{\link{ppm}()}, not fixed in
  \code{MultiStraussHard()}. 
}
\seealso{
  \code{\link{ppm}},
  \code{\link{pairwise.family}},
  \code{\link{ppm.object}},
  \code{\link{MultiStrauss}},
  \code{\link{MultiHard}},
  \code{\link{Strauss}}
}
\examples{
   r <- matrix(3, nrow=2,ncol=2)
   h <- matrix(c(1,2,2,1), nrow=2,ncol=2)
   MultiStraussHard(r,h)
   # prints a sensible description of itself
   r <- 0.04 * matrix(c(1,2,2,1), nrow=2,ncol=2)
   h <- 0.02 * matrix(c(1,NA,NA,1), nrow=2,ncol=2)
   X <- amacrine
   \testonly{
       X <- X[owin(c(0,0.8), c(0,1))]
   }
   fit <- ppm(X ~1, MultiStraussHard(r,h))
   # fit stationary multitype hardcore Strauss process to `amacrine'
}
\section{Warnings}{
  In order that \code{\link{ppm}} can fit the multitype/hard core
  Strauss model correctly to a point pattern \code{X}, this pattern
  must be marked, with \code{markformat} equal to \code{vector}
  and the mark vector \code{marks(X)} must be a factor.  If the
  argument \code{types} is specified it is interpreted as a set of
  factor levels and this set must equal \code{levels(marks(X))}.
}
\section{Changed Syntax}{
  Before \pkg{spatstat} version \code{1.37-0},
  the syntax of this function was different:
  \code{MultiStraussHard(types=NULL, iradii, hradii)}.
  The new code attempts to handle the old syntax as well.
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  ,
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{models}
