\name{hybrid.family}
\alias{hybrid.family}
\title{
  Hybrid Interaction Family
}
\description{
  An object describing the family of all hybrid interactions.
}
\details{
  \bold{Advanced Use Only!}
  
  This structure would not normally be touched by
  the user. It describes 
  the family of all hybrid point process models.
 
  If you need to create a specific hybrid interaction model for use in 
  modelling, use the function \code{\link{Hybrid}}.
 
  Anyway, \code{hybrid.family} is an object of class \code{"isf"}
  containing a function \code{hybrid.family$eval} for
  evaluating the sufficient statistics of any hybrid interaction
  point process model.
}
\seealso{
  Use \code{\link{Hybrid}} to make hybrid interactions.
  
  Other families: 
  \code{\link{pairwise.family}},
  \code{\link{pairsat.family}},
  \code{\link{ord.family}},
  \code{\link{inforder.family}}.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
