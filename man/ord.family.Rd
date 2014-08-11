\name{ord.family}
\alias{ord.family}
\title{Ord Interaction Process Family}
\description{
  An object describing the family of all Ord interaction point processes
}
\details{
  \bold{Advanced Use Only!}
  
  This structure would not normally be touched by
  the user. It describes 
  the family of point process models introduced by Ord (1977).
 
  If you need to create a specific Ord-type model for use in 
  analysis, use the function \code{\link{OrdThresh}} or
  \code{\link{Ord}}.
 
  Anyway, \code{ord.family} is an object of class \code{"isf"}
  containing a function \code{ord.family$eval} for
  evaluating the sufficient statistics of any Ord type
  point process model taking an exponential family form. 
} 
\seealso{
  \code{\link{pairwise.family}},
  \code{\link{pairsat.family}},
  \code{\link{Poisson}},
  \code{\link{Pairwise}},
  \code{\link{PairPiece}},
  \code{\link{Strauss}},
  \code{\link{StraussHard}},
  \code{\link{Softcore}},
  \code{\link{Geyer}},
  \code{\link{SatPiece}},
  \code{\link{Saturated}},
  \code{\link{Ord}},
  \code{\link{OrdThresh}}
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
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
