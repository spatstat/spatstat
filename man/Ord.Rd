\name{Ord}
\alias{Ord}
\title{Generic Ord Interaction model}
\description{
Creates an instance of an Ord-type interaction point process model
which can then be fitted to point pattern data.
}
\usage{
  Ord(pot, name)
}
\arguments{
  \item{pot}{An S language function giving the user-supplied
    interaction potential.}
  \item{name}{Character string.}
}
\value{
  An object of class \code{"interact"}
  describing the interpoint interaction
  structure of a point process. 
}
\details{
  Ord's point process model (Ord, 1977) is a Gibbs point process
  of infinite order. Each point \eqn{x_i}{x[i]} in the point pattern
  \eqn{x} contributes a factor \eqn{g(a_i)}{g(a[i])} where
  \eqn{a_i = a(x_i, x)}{a[i] = a(x[i], x)} is the area of the
  tile associated with \eqn{x_i}{x[i]}
  in the Dirichlet tessellation of \eqn{x}.

  Ord (1977) proposed fitting this model to forestry data
  when \eqn{g(a)} has a simple ``threshold'' form. That model is
  implemented in our function \code{\link{OrdThresh}}.
  The present function \code{Ord} implements the case of a
  completely general Ord potential \eqn{g(a)}
  specified as an S language function \code{pot}.

  This is experimental. 
}
\references{
  Baddeley, A. and Turner, R. (2000)
  Practical maximum pseudolikelihood for spatial point patterns.
  \emph{Australian and New Zealand Journal of Statistics}
  \bold{42}, 283--322.

  Ord, J.K. (1977) 
  Contribution to the discussion of Ripley (1977).

  Ord, J.K. (1978) 
  How many trees in a forest?
  \emph{Mathematical Scientist} \bold{3}, 23--33.

  Ripley, B.D. (1977)
  Modelling spatial patterns (with discussion).
  \emph{Journal of the Royal Statistical Society, Series B},
  \bold{39}, 172 -- 212.

}
\seealso{
  \code{\link{ppm}},
  \code{\link{ppm.object}},
  \code{\link{OrdThresh}}
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{models}

