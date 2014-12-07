\name{hierpair.family}
\alias{hierpair.family}
\title{Hierarchical Pairwise Interaction Process Family}
\description{
  An object describing the family of all hierarchical pairwise interaction Gibbs
  point processes.
}
\details{
  \bold{Advanced Use Only!}
  
  This structure would not normally be touched by
  the user. It describes 
  the hierarchical pairwise interaction family of point process models.

  Anyway, \code{hierpair.family} is an object of class \code{"isf"}
  containing a function \code{hierpair.family$eval} for
  evaluating the sufficient statistics of any hierarchical pairwise interaction
  point process model taking an exponential family form. 
} 
\seealso{
  Other families: 
  \code{\link[spatstat]{pairwise.family}},
  \code{\link[spatstat]{pairsat.family}},
  \code{\link[spatstat]{ord.family}},
  \code{\link[spatstat]{inforder.family}}.

  Hierarchical Strauss interaction:
  \code{\link{HierStrauss}}.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak.
}
\keyword{spatial}
\keyword{models}
