\name{inforder.family}
\alias{inforder.family}
\title{Infinite Order Interaction Family}
\description{
  An object describing the family of all Gibbs point processes
  with infinite interaction order.
}
\details{
  \bold{Advanced Use Only!}
  
  This structure would not normally be touched by
  the user. It describes the interaction structure
  of Gibbs point processes which have infinite order of interaction,
  such as the area-interaction process \cite{\link{AreaInter}}.
 
  Anyway, \code{inforder.family} is an object of class \code{"isf"}
  containing a function \code{inforder.family$eval} for
  evaluating the sufficient statistics of a Gibbs
  point process model taking an exponential family form. 
} 
\seealso{
  \code{\link{AreaInter}} to create the area interaction process
  structure.
  
  Other families:
  \code{\link{pairwise.family}},
  \code{\link{pairsat.family}},
  \code{\link{ord.family}}.

  
}
\references{
  Baddeley, A. and Turner, R. (2000)
  Practical maximum pseudolikelihood for spatial point patterns.
  \emph{Australian and New Zealand Journal of Statistics}
  \bold{42}, 283--322.
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{models}
