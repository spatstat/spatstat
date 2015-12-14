\name{branchlabelfun}
\alias{branchlabelfun}
\title{
  Tree Branch Membership Labelling Function
}
\description{
  Creates a function which returns the tree branch membership label
  for any location on a linear network.
}
\usage{
  branchlabelfun(L, root = 1)
}
\arguments{
  \item{L}{
    Linear network (object of class \code{"linnet"}).
    The network must have no loops.
  }
  \item{root}{
    Root of the tree. An integer index identifying
    which point in \code{vertices(L)} is the root of the tree.
  }
}
\details{
  The linear network \code{L}
  must be an acyclic graph (i.e. must not contain any loops) so that it
  can be interpreted as a tree. 
  
  The result of \code{f <- branchlabelfun(L, root)} is
  a function \code{f} which gives,
  for each location on the linear network \code{L},
  the tree branch label at that location.

  Tree branch labels are explained in \code{\link{treebranchlabels}}.
  
  The result \code{f} also belongs to the class \code{"linfun"}.
  It can be called using several different kinds of data,
  as explained in the help for \code{\link{linfun}}.
  The values of the function are character strings.
}
\value{
  A function (of class \code{"linfun"}).
}
\author{
Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}


Rolf Turner \email{r.turner@auckland.ac.nz}

and Ege Rubak \email{rubak@math.aau.dk}

}
\seealso{
  \code{\link{treebranchlabels}},
  \code{\link{linfun}}
}
\examples{
  # make a simple tree
  m <- simplenet$m
  m[8,10] <- m[10,8] <- FALSE
  L <- linnet(vertices(simplenet), m)
  # make function
  f <- branchlabelfun(L, 1)
  plot(f)
  X <- runiflpp(5, L)
  f(X)
}
\keyword{spatial}
\keyword{math}

