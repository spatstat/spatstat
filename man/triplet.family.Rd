\name{triplet.family}
\alias{triplet.family}
\title{Triplet Interaction Family}
\description{
  An object describing the family of all Gibbs point processes
  with interaction order equal to 3.
}
\details{
  \bold{Advanced Use Only!}
  
  This structure would not normally be touched by
  the user. It describes the interaction structure
  of Gibbs point processes which have infinite order of interaction,
  such as the triplet interaction process \cite{\link{Triplets}}.
 
  Anyway, \code{triplet.family} is an object of class \code{"isf"}
  containing a function \code{triplet.family$eval} for
  evaluating the sufficient statistics of a Gibbs
  point process model taking an exponential family form. 
} 
\seealso{
  \code{\link{Triplets}} to create the triplet interaction process
  structure.
  
  Other families:
  \code{\link{pairwise.family}},
  \code{\link{pairsat.family}},
  \code{\link{inforder.family}},
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
