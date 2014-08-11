\name{intensity}
\alias{intensity}
\title{
  Intensity of a Dataset or a Model
}
\description{
  Generic function for computing the intensity of a spatial dataset
  or spatial point process model.
}
\usage{
intensity(X, ...)
}
\arguments{
  \item{X}{
    A spatial dataset or a spatial point process model.
  }
  \item{\dots}{
    Further arguments depending on the class of \code{X}.
  }
}
\details{
  This is a generic function for computing the intensity of a spatial dataset
  or spatial point process model. There are methods for point patterns
  (objects of class \code{"ppp"}) and fitted point process models
  (objects of class \code{"ppm"}).
  
  The empirical intensity of a dataset is the average density
  (the average amount of \sQuote{stuff} per unit area or volume).
  The empirical intensity of a point pattern is computed by the
  method \code{\link{intensity.ppp}}.

  The theoretical intensity of a stochastic model is the expected density
  (expected  amount of \sQuote{stuff} per unit area or volume).
  The theoretical intensity of a fitted point process model is computed by the
  method \code{\link{intensity.ppm}}.
}
\value{
  Usually a numeric value or vector.
}
\seealso{
  \code{\link{intensity.ppp}},
  \code{\link{intensity.ppm}}.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
