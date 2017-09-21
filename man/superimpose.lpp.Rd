\name{superimpose.lpp}  
\alias{superimpose.lpp}
\title{Superimpose Several Point Patterns on Linear Network}
\description{
  Superimpose any number of point patterns on the same linear network.
}
\usage{
  \method{superimpose}{lpp}(\dots, L=NULL)
}
\arguments{
  \item{\dots}{
    Any number of arguments, each of which represents a point
    pattern on the same linear network.
    Each argument can be either an object of class \code{"lpp"},
    giving both the spatial coordinates of the points and the
    linear network, or a \code{list(x,y)} or \code{list(x,y,seg,tp)}
    giving just the spatial coordinates of the points.
  }
  \item{L}{
    Optional. The linear network.
    An object of class \code{"linnet"}.
    This argument is required if none of the other arguments
    is of class \code{"lpp"}. 
  }
}
\value{
  An object of class \code{"lpp"} representing the combined point
  pattern on the linear network.
}
\details{
  This function is used to superimpose several point patterns
  on the same linear network. It is a method for the generic
  function \code{\link{superimpose}}.

  Each of the arguments \code{\dots}
  can be either a point pattern on a linear network
  (object of class \code{"lpp"} 
  giving both the spatial coordinates of the points and the
  linear network), or a \code{list(x,y)} or \code{list(x,y,seg,tp)}
  giving just the spatial coordinates of the points.
  These arguments must represent point patterns on the \emph{same}
  linear network.

  The argument \code{L} is an alternative way to specify the linear
  network, and is required if none of the arguments \code{\dots} is an
  object of class \code{"lpp"}.

  The arguments \code{\dots} may be \emph{marked} patterns.
  The marks of each component pattern must have the same format.
  Numeric and character marks may be ``mixed''.  If there is such
  mixing then the numeric marks are coerced to character in the
  combining process. If the mark structures are all data frames,
  then these data frames must have the same number of columns and
  identical column names.

  If the arguments \code{\dots} are given in the form \code{name=value},
  then the \code{name}s will be used as an extra column of marks
  attached to the elements of the corresponding patterns.
}
\seealso{
  \code{\link{superimpose}}
}
\examples{
  X <- rpoislpp(5, simplenet)
  Y <- rpoislpp(10, simplenet)
  superimpose(X,Y) # not marked
  superimpose(A=X, B=Y) # multitype with types A and B
}
\author{\adrian
  
  
  \rolf
  
  \ege
  
  and Greg McSwiggan.
}
\keyword{spatial}
\keyword{manip}
