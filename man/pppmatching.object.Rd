\name{pppmatching.object}
\alias{pppmatching.object} %DoNotExport
\title{Class of Point Matchings}
\description{
  A class \code{"pppmatching"} to represent a matching of two planar
  point patterns.
  Optionally includes information about the construction of the matching
  and its associated distance between the point patterns.
}
\details{
  This class represents a (possibly weighted and incomplete) matching
  between two planar point patterns (objects of class \code{"ppp"}).

  A matching can be thought of as a bipartite weighted graph where
  the vertices are given by the two point patterns and edges of positive
  weights are drawn each time a point of the first point pattern is
  "matched" with a point of the second point pattern.

  If \code{m} is an object of type \code{pppmatching}, it contains the
  following elements
  \tabular{ll}{
    \code{pp1, pp2} \tab the two point patterns to be matched (vertices) \cr
    \code{matrix} \tab a matrix specifying which points are matched \cr
                  \tab and with what weights (edges) \cr
    \code{type} \tab (optional) a character string for the type of \cr
                 \tab the matching (one of \code{"spa"}, \code{"ace"} or \code{"mat"}) \cr
    \code{cutoff} \tab (optional) cutoff value for interpoint distances \cr
    \code{q} \tab (optional) the order for taking averages of \cr
             \tab interpoint distances \cr
    \code{distance} \tab (optional) the distance associated with the matching
  }

  The element \code{matrix} is a "generalized adjacency matrix".
  The numbers of rows
  and columns match the cardinalities of the first and second point
  patterns, respectively. The \code{[i,j]}-th entry is positive if
  the \code{i}-th point of \code{X} and the \code{j}-th point of
  \code{Y} are matched (zero otherwise) and its value then gives
  the corresponding weight of the match. For an unweighted matching
  all the weights are set to \eqn{1}.

  The optional elements are for saving details about matchings in the context of
  optimal point matching techniques. \code{type} can be one of \code{"spa"} (for
  "subpattern assignment"), \code{"ace"} (for "assignment only if
  cardinalities differ") or \code{"mat"} (for "mass transfer"). \code{cutoff}
  is a positive numerical value that specifies the maximal interpoint distance and
  \code{q} is a value in \eqn{[1,\infty]}{[1,Inf]} that gives the order of the average
  applied to the interpoint distances. See the help files for \code{\link{pppdist}}
  and \code{\link{matchingdist}} for detailed information about these elements.
  
  Objects of class \code{"pppmatching"} may be created by the function 
  \code{\link{pppmatching}}, and are most commonly obtained as output of the
  function \code{\link{pppdist}}. There are methods \code{plot}, \code{print} and
  \code{summary} for this class.
}
\author{
  Dominic Schuhmacher
  \email{dominic.schuhmacher@stat.unibe.ch}
  \url{http://www.dominic.schuhmacher.name}
}
\seealso{
  \code{\link{matchingdist}}
  \code{\link{pppmatching}}
}
\examples{
  # a random complete unweighted matching
  X <- runifpoint(10)
  Y <- runifpoint(10)
  am <- r2dtable(1, rep(1,10), rep(1,10))[[1]]
        # generates a random permutation matrix
  m <- pppmatching(X, Y, am)
  summary(m)
  m$matrix
  \dontrun{
    plot(m)
  }

  # an optimal complete unweighted matching
  m2 <- pppdist(X,Y)
  summary(m2)
  m2$matrix
  \dontrun{
    plot(m2)
  }
}
\keyword{spatial}
\keyword{attribute}
