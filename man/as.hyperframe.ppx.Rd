\name{as.hyperframe.ppx}
\Rdversion{1.1}
\alias{as.hyperframe.ppx}
\alias{as.data.frame.ppx}
\alias{as.matrix.ppx}
\title{
  Extract coordinates and marks of multidimensional point pattern
}
\description{
  Given any kind of spatial or space-time point pattern,
  extract the coordinates and marks of the points.
}
\usage{
\method{as.hyperframe}{ppx}(x, ...)
\method{as.data.frame}{ppx}(x, ...)
\method{as.matrix}{ppx}(x, ...)
}
\arguments{
  \item{x}{
    A general multidimensional space-time point pattern 
    (object of class \code{"ppx"}).
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  An object of class \code{"ppx"} (see \code{\link{ppx}})
  represents a marked point pattern
  in multidimensional space and/or time. There may be any
  number of spatial coordinates, any number of temporal coordinates,
  and any number of mark variables. The individual marks may be
  atomic (numeric values, factor values, etc) or objects of any kind.

  The function \code{as.hyperframe.ppx} extracts the coordinates
  and the marks as a \code{"hyperframe"} (see
  \code{\link{hyperframe}}) with one row of data for each
  point in the pattern. This is a method for the generic
  function \code{\link{as.hyperframe}}.

  The function \code{as.data.frame.ppx} discards those mark variables
  which are not atomic values, and extracts the coordinates
  and the remaining marks as a \code{data.frame}
  with one row of data for each
  point in the pattern. This is a method for the generic
  function \code{\link{as.data.frame}}.

  Finally \code{as.matrix(x)} is equivalent to
  \code{as.matrix(as.data.frame(x))} for an object of class
  \code{"ppx"}. Be warned that, if there are any columns of non-numeric
  data (i.e. if there are mark variables that are factors),
  the result will be a matrix of character values.
}
\value{
  A \code{hyperframe}, \code{data.frame} or \code{matrix} as appropriate.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{ppx}},
  \code{\link{hyperframe}},
  \code{\link{as.hyperframe}}.
}
\examples{
   df <- data.frame(x=runif(4),y=runif(4),t=runif(4))
   X <- ppx(data=df, coord.type=c("s","s","t"))
   as.data.frame(X)
   val <- runif(4)
   E <- lapply(val, function(s) { rpoispp(s) })
   hf <- hyperframe(t=val, e=as.listof(E))
   Z <- ppx(data=hf, domain=c(0,1))
   as.hyperframe(Z)
   as.data.frame(Z)
}
\keyword{spatial}
\keyword{manip}
