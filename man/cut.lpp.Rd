\name{cut.lpp}
\alias{cut.lpp}
\title{Classify Points in a Point Pattern on a Network}
\description{
  For a point pattern on a linear network, 
  classify the points into distinct types
  according to the numerical marks in the pattern, or according to
  another variable.
}
\usage{
  \method{cut}{lpp}(x, z=marks(x), ...)
}
\arguments{
  \item{x}{
    A point pattern on a linear network
    (object of class \code{"lpp"}).
  }
  \item{z}{
    Data determining the classification. A numeric vector,
    a factor, a pixel image on a linear network (class \code{"linim"}),
    a function on a linear network (class \code{"linfun"}),
    a tessellation on a linear network (class \code{"lintess"}), a string
    giving the name of a column of marks, or one of the coordinate
    names \code{"x"}, \code{"y"}, \code{"seg"} or \code{"tp"}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{cut.default}}.
    They determine the breakpoints for the mapping from numerical values
    in \code{z} to factor values in the output.
    See \code{\link{cut.default}}.
  }
} 
\value{
  A multitype point pattern on the same linear network,
  that is, a point pattern object
  (of class \code{"lpp"}) with a \code{marks} vector that is a factor.
}
\details{
  This function has the effect of classifying each point in the point
  pattern \code{x} into one of several possible types. The
  classification is based on the dataset \code{z}, which may be either
  \itemize{
    \item
    a factor (of length equal to the number of points in \code{z})
    determining the classification of each point in \code{x}.
    Levels of the factor determine the classification.
    \item
    a numeric vector (of length equal to the number of points in
    \code{z}). The range of values of \code{z} will be divided into
    bands (the number of bands is determined by \code{\dots})
    and \code{z} will be converted to a factor using
    \code{\link{cut.default}}.
    \item
    a pixel image on a network (object of class \code{"linim"}).
    The value of \code{z} at each point of \code{x} will be
    used as the classifying variable.
    \item
    a function on a network (object of class \code{"linfun"}, see
    \code{\link{linfun}}). 
    The value of \code{z} at each point of \code{x} will be
    used as the classifying variable.
    \item
    a tessellation on a network (object of class \code{"lintess"}, see
    \code{\link{lintess}}). Each point of \code{x} will be classified
    according to the tile of the tessellation into which it falls.
    \item
    a character string, giving the name of one of the columns
    of \code{marks(x)}, if this is a data frame.
    \item
    a character string identifying one of the coordinates:
    the spatial coordinates
    \code{"x"}, \code{"y"} or the segment identifier \code{"seg"}
    or the fractional coordinate along the segment, \code{"tp"}.
  }
  The default is to take \code{z} to be the vector of marks in
  \code{x} (or the first column in the data frame of marks of \code{x},
  if it is a data frame). If the marks are numeric, then the range of values
  of the numerical marks is divided into several intervals, and each
  interval is associated with a level of a factor. 
  The result is a
  marked point pattern, on the same linear network,
  with the same point locations as
  \code{x}, but with the numeric mark of each point discretised
  by replacing it by the factor level.
  This is a convenient way to transform a marked point pattern
  which has numeric marks into a multitype point pattern,
  for example to plot it or analyse it. See the examples.

  To select some points from \code{x}, use the subset operators
  \code{\link{[.lpp}} or \code{\link{subset.lpp}} instead.
}
\seealso{
  \code{\link{cut}},
  \code{\link{lpp}},
  \code{\link{lintess}},
  \code{\link{linfun}},
  \code{\link{linim}}
}
\examples{
  X <- runiflpp(20, simplenet)
  f <- linfun(function(x,y,seg,tp) { x }, simplenet)
  plot(cut(X, f, breaks=4))
  plot(cut(X, "x", breaks=4))
  plot(cut(X, "seg"))
}

\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{methods}
