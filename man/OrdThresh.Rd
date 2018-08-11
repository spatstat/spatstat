\name{OrdThresh}
\alias{OrdThresh}
\title{Ord's Interaction model}
\description{
Creates an instance of Ord's point process model
which can then be fitted to point pattern data.
}
\usage{
  OrdThresh(r)
}
\arguments{
  \item{r}{Positive number giving the threshold value
    for Ord's model.}
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
  in the Dirichlet tessellation of \eqn{x}. The function \eqn{g} is
  simply \eqn{g(a) = 1} if \eqn{a \ge r}{a >= r} and
  \eqn{g(a) = \gamma < 1}{g(a) = gamma < 1} if \eqn{a < r}{a < r},
  where \eqn{r} is called the threshold value.

  This function creates an instance of Ord's model with a given
  value of \eqn{r}. It can then be fitted to point process data
  using \code{\link{ppm}}. 
}
\seealso{
  \code{\link{ppm}},
  \code{\link{ppm.object}}
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
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{models}

