\name{markstat}
\alias{markstat}
\title{Summarise Marks in Every Neighbourhood in a Point Pattern}
\description{
  Visit each point in a point pattern, find the neighbouring points,
  and summarise their marks
}
\usage{
   markstat(X, fun, N=NULL, R=NULL, \dots)
}
\arguments{
  \item{X}{
    A marked point pattern. 
    An object of class \code{"ppp"}.
  }
  \item{fun}{
    Function to be applied to the vector of marks.
  }
  \item{N}{
    Integer. If this argument is present,
    the neighbourhood of a point of \code{X} is defined to consist of the
    \code{N} points of \code{X} which are closest to it.
  }
  \item{R}{
    Nonnegative numeric value. If this argument is present,
    the neighbourhood of a point of \code{X} is defined to consist of
    all points of \code{X} which lie within a distance \code{R}
    of it.
  }
  \item{\dots}{
    extra arguments passed to the function \code{fun}.
    They must be given in the form \code{name=value}.
  }
}
\value{
  Similar to the result of \code{\link{apply}}.
  if each call to \code{fun} returns a single numeric value,
  the result is a vector of dimension \code{npoints(X)}, the number of points
  in \code{X}.
  If each call to \code{fun} returns a vector of the same length
  \code{m}, then the result is a matrix of dimensions \code{c(m,n)};
  note the transposition of the indices, as usual for the family of
  \code{apply} functions.
  If the calls to \code{fun} return vectors of different lengths,
  the result is a list of length \code{npoints(X)}. 
}
\details{
  This algorithm visits each point in the point pattern \code{X},
  determines which points of \code{X} are ``neighbours'' of the current
  point, extracts the marks of these neighbouring points,
  applies the function \code{fun} to the marks,
  and collects the value or values returned by \code{fun}.

  The definition of ``neighbours'' depends on the arguments
  \code{N} and \code{R}, exactly one of which must be given.

  If \code{N} is given, then the neighbours of the current
  point are the \code{N} points of \code{X} which are closest to
  the current point (including the current point itself).
  If \code{R} is given, then the neighbourhood of the current point
  consists of all points of \code{X} which lie closer than a distance \code{R}
  from the current point.

  Each point of \code{X} is visited; the neighbourhood
  of the current point is determined; the marks of these points
  are extracted as a vector \code{v}; then the function 
  \code{fun} is called as:
  
  \code{fun(v, \dots)}

  where \code{\dots} are the arguments passed from the call to
  \code{markstat}.

  The results of each call to \code{fun} are collected and returned
  according to the usual rules for \code{\link{apply}} and its
  relatives. See the section on \bold{Value}.

  This function is just a convenient wrapper for a common use of the
  function \code{\link{applynbd}}. For more complex tasks,
  use \code{\link{applynbd}}. 
  To simply tabulate the marks in every \code{R}-neighbourhood, use
  \code{\link{marktable}}.
}
\seealso{
  \code{\link{applynbd}},
  \code{\link{marktable}},
  \code{\link{ppp.object}},
  \code{\link{apply}}
}
\examples{
  trees <- longleaf
  \testonly{
     trees <- trees[seq(1, npoints(trees), by=6)]
  }

  # average diameter of 5 closest neighbours of each tree
  md <- markstat(trees, mean, N=5)

  # range of diameters of trees within 10 metre radius
  rd <- markstat(trees, range, R=10)
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au},
 Rolf Turner \email{r.turner@auckland.ac.nz}
 and Ege Rubak \email{rubak@math.aau.dk}.
}
\keyword{spatial}
\keyword{programming}
 
