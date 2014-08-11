\name{applynbd}
\alias{applynbd}
\title{Apply Function to Every Neighbourhood in a Point Pattern}
\description{
  Visit each point in a point pattern, find the neighbouring points,
  and apply a given function to them.
}
\usage{
   applynbd(X, FUN, N, R, criterion, exclude=FALSE, \dots)
}
\arguments{
  \item{X}{
    Point pattern. 
    An object of class \code{"ppp"},
    or data which can be converted into 
    this format by \code{\link{as.ppp}}.
  }
  \item{FUN}{
    Function to be applied to each neighbourhood.
    The arguments of \code{FUN} are described under \bold{Details}.
  }
  \item{N}{
    Integer. If this argument is present,
    the neighbourhood of a point of \code{X} is defined to consist of the
    \code{N} points of \code{X} which are closest to it.
    This argument is incompatible with \code{R} and \code{criterion}.
  }
  \item{R}{
    Nonnegative numeric value. If this argument is present,
    the neighbourhood of a point of \code{X} is defined to consist of
    all points of \code{X} which lie within a distance \code{R}
    of it.
    This argument is incompatible with \code{N} and \code{criterion}.
  }
  \item{criterion}{
    Function. If this argument is present, 
    the neighbourhood of a point of \code{X} is determined by
    evaluating this function. See under \bold{Details}.
    This argument is incompatible with \code{N} and \code{R}.
  }
  \item{exclude}{
    Logical. If \code{TRUE} then the point currently being visited
    is excluded from its own neighbourhood.
  }
  \item{\dots}{
    extra arguments passed to the function \code{FUN}.
    They must be given in the form \code{name=value}.
  }
}
\value{
  Similar to the result of \code{\link{apply}}.
  If each call to \code{FUN} returns a single numeric value,
  the result is a vector of dimension \code{X$n}, the number of points
  in \code{X}.
  If each call to \code{FUN} returns a vector of the same length
  \code{m}, then the result is a matrix of dimensions \code{c(m,n)};
  note the transposition of the indices, as usual for the family of
  \code{apply} functions.
  If the calls to \code{FUN} return vectors of different lengths,
  the result is a list of length \code{X$n}. 
}
\details{
  This is an analogue of \code{\link{apply}}
  for point patterns. It visits each point in the point pattern \code{X},
  determines which points of \code{X} are ``neighbours'' of the current
  point, applies the function \code{FUN} to this neighbourhood,
  and collects the values returned by \code{FUN}.

  The definition of ``neighbours'' depends on the arguments
  \code{N}, \code{R} and \code{criterion}, exactly one of which
  must be given.
  Also the argument \code{exclude} determines whether
  the current point is excluded from its own neighbourhood.

  \itemize{
    \item 
    If \code{N} is given, then the neighbours of the current
    point are the \code{N} points of \code{X} which are closest to
    the current point (including the current point itself
    unless \code{exclude=TRUE}).
    \item
    If \code{R} is given, then the neighbourhood of the current point
    consists of all points of \code{X} which lie closer than a distance \code{R}
    from the current point.
    \item 
    If \code{criterion} is given, then it must be a function
    with two arguments \code{dist} and \code{drank} which will be
    vectors of equal length.
    The interpretation is that \code{dist[i]} will be the
    distance of a point from the current point, and
    \code{drank[i]} will be the rank of that distance (the three points
    closest to the current point will have rank 1, 2 and 3).
    This function must return a logical vector of the same length
    as \code{dist} and \code{drank} whose \code{i}-th entry is
    \code{TRUE} if the corresponding point should be included in
    the neighbourhood. See the examples below.
  }

  When \code{applynbd} is executed, 
  each point of \code{X} is visited, and the following happens
  for each point:
  \itemize{
    \item
    the neighbourhood of the current point is determined according
    to the chosen rule, and stored as a point pattern \code{Y};
    \item
    the function \code{FUN} is called as:
  
    \code{FUN(Y=Y, current=current, dists=dists, dranks=dranks, \dots)}

    where \code{current} is the location of the current point
    (in a format explained below),
    \code{dists} is a vector of distances from the current
    point to each of the points in \code{Y}, 
    \code{dranks} is a vector of the ranks of these distances
    with respect to the full point pattern \code{X},
    and \code{\dots} are the arguments passed from the call to
    \code{applynbd};
    \item
    The result of the call to \code{FUN} is stored.
  }
  The results of each call to \code{FUN} are collected and returned
  according to the usual rules for \code{\link{apply}} and its
  relatives. See the \bold{Value} section of this help file.

  The format of the argument \code{current} is as follows.
  If \code{X} is an unmarked point pattern, then \code{current} is a
  vector of length 2 containing the coordinates of the current point.
  If \code{X} is marked, then \code{current} is a point pattern
  containing exactly one point, so that \code{current$x} is its
  \eqn{x}-coordinate and \code{current$marks} is its mark value.
  In either case, the coordinates of the current point can be referred to as
  \code{current$x} and \code{current$y}.

  Note that \code{FUN} will be called exactly as described above,
  with each argument named explicitly. Care is required when writing the
  function \code{FUN} to ensure that
  the arguments will match up. See the Examples.

  See \code{\link{markstat}} for a common use of this function.

  To simply tabulate the marks in every \code{R}-neighbourhood, use
  \code{\link{marktable}}.
}
\seealso{
  \code{\link{ppp.object}},
  \code{\link{apply}},
  \code{\link{markstat}},
  \code{\link{marktable}}
}
\examples{
  data(redwood)

  # count the number of points within radius 0.2 of each point of X
  nneighbours <- applynbd(redwood, R=0.2, function(Y, ...){Y$n-1})
  # equivalent to:
  nneighbours <- applynbd(redwood, R=0.2, function(Y, ...){Y$n}, exclude=TRUE)

  # compute the distance to the second nearest neighbour of each point
  secondnndist <- applynbd(redwood, N = 2,
                           function(dists, ...){max(dists)},
                           exclude=TRUE)

  # marked point pattern
  data(longleaf)
  \testonly{
	# smaller dataset
	longleaf <- longleaf[seq(1, longleaf$n, by=80)]
  }
  # compute the median of the marks of all neighbours of a point
  # (see also 'markstat')
  dbh.med <- applynbd(longleaf, R=90, exclude=TRUE,
                 function(Y, ...) { median(Y$marks)})


  # ANIMATION explaining the definition of the K function
  # (arguments `fullpicture' and 'rad' are passed to FUN)

  \dontrun{
  showoffK <- function(Y, current, dists, dranks, fullpicture,rad) { 
	plot(fullpicture, main="")
	points(Y, cex=2)
        u <- current
	points(u$x,u$y,pch="+",cex=3)
	theta <- seq(0,2*pi,length=100)
	polygon(u$x+ rad * cos(theta),u$y+rad*sin(theta))
	text(u$x+rad/3,u$y+rad/2,Y$n,cex=3)
	Sys.sleep(if(runif(1) < 0.1) 1.5 else 0.3)
	return(Y$n - 1)
  }
  applynbd(redwood, R=0.2, showoffK, fullpicture=redwood, rad=0.2, exclude=TRUE)

  # animation explaining the definition of the G function

  showoffG <- function(Y, current, dists, dranks, fullpicture) { 
	plot(fullpicture, main="")
	points(Y, cex=2)
        u <- current
	points(u[1],u[2],pch="+",cex=3)
	v <- c(Y$x[1],Y$y[1])
	segments(u[1],u[2],v[1],v[2],lwd=2)
	w <- (u + v)/2
	nnd <- dists[1]
	text(w[1],w[2],round(nnd,3),cex=2)
	Sys.sleep(if(runif(1) < 0.1) 1.5 else 0.3)
	return(nnd)
  }

  data(cells)
  applynbd(cells, N=1, showoffG, exclude=TRUE, fullpicture=cells)
  }
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{programming}
\keyword{iteration}
 
