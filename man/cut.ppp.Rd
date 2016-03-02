\name{cut.ppp}
\alias{cut.ppp}
\title{Classify Points in a Point Pattern}
\description{
  Classifies the points in a point pattern into distinct types
  according to the numerical marks in the pattern, or according to
  another variable.
}
\usage{
  \method{cut}{ppp}(x, z=marks(x), ...)
}
\arguments{
  \item{x}{
    A two-dimensional point pattern.
    An object of class \code{"ppp"}.
  }
  \item{z}{
    Data determining the classification. A numeric vector,
    a factor, a pixel image, a window, a tessellation, or a string
    giving the name of a column of marks.
  }
  \item{\dots}{
    Arguments passed to \code{\link{cut.default}}.
    They determine the breakpoints for the mapping from numerical values
    in \code{z} to factor values in the output.
    See \code{\link{cut.default}}.
  }
} 
\value{
  A multitype point pattern, that is, a point pattern object
  (of class \code{"ppp"}) with a \code{marks} vector that is a factor.
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
    a pixel image (object of class \code{"im"}).
    The value of \code{z} at each point of \code{x} will be
    used as the classifying variable.
    \item
    a tessellation (object of class \code{"tess"}, see
    \code{\link{tess}}). Each point of \code{x} will be classified
    according to the tile of the tessellation into which it falls.
    \item
    a window (object of class \code{"owin"}).
    Each point of \code{x} will be classified
    according to whether it falls inside or outside this window.
    \item
    a character string, giving the name of one of the columns
    of \code{marks(x)}, if this is a data frame.
  }
  The default is to take \code{z} to be the vector of marks in
  \code{x} (or the first column in the data frame of marks of \code{x},
  if it is a data frame). If the marks are numeric, then the range of values
  of the numerical marks is divided into several intervals, and each
  interval is associated with a level of a factor. 
  The result is a
  marked point pattern, with the same window and point locations as
  \code{x}, but with the numeric mark of each point discretised
  by replacing it by the factor level.
  This is a convenient way to transform a marked point pattern
  which has numeric marks into a multitype point pattern,
  for example to plot it or analyse it. See the examples.

  To select some points from a point pattern, use the subset operator
  \code{\link{[.ppp}} instead.
}
\seealso{
  \code{\link{cut}},
  \code{\link{ppp.object}},
  \code{\link{tess}}
}
\examples{
 # (1) cutting based on numeric marks of point pattern
 
 trees <- longleaf
 # Longleaf Pines data
 # the marks are positive real numbers indicating tree diameters.

 \testonly{
	# smaller dataset
	trees <- trees[seq(1, npoints(trees), by=80)]
 }
 \dontrun{
 plot(trees)
 }

 # cut the range of tree diameters into three intervals
 long3 <- cut(trees, breaks=3)
 \dontrun{
 plot(long3)
 }

 # adult trees defined to have diameter at least 30 cm
 long2 <- cut(trees, breaks=c(0,30,100), labels=c("Sapling", "Adult"))
 plot(long2)
 plot(long2, cols=c("green","blue"))

 # (2) cutting based on another numeric vector
 # Divide Swedish Pines data into 3 classes
 # according to nearest neighbour distance

 swedishpines
 plot(cut(swedishpines, nndist(swedishpines), breaks=3))

 # (3) cutting based on tessellation
 # Divide Swedish Pines study region into a 4 x 4 grid of rectangles
 # and classify points accordingly

 tes <- tess(xgrid=seq(0,96,length=5),ygrid=seq(0,100,length=5))
 plot(cut(swedishpines, tes))
 plot(tes, lty=2, add=TRUE)

 # (4) multivariate marks
 finpines
 cut(finpines, "height", breaks=4)
}

\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au},
 Rolf Turner \email{r.turner@auckland.ac.nz}
 and Ege Rubak \email{rubak@math.aau.dk}.
}
\keyword{spatial}
\keyword{methods}
