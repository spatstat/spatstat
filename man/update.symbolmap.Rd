\name{update.symbolmap}
\alias{update.symbolmap}
\title{
  Update a Graphics Symbol Map.
}
\description{
  This command updates the \code{object} using the
  arguments given.
}
\usage{
\method{update}{symbolmap}(object, \dots)
}
\arguments{
  \item{object}{
    Graphics symbol map (object of class \code{"symbolmap"}).
  }
  \item{\dots}{
    Additional or replacement arguments to \code{\link{symbolmap}}.
  }
}
\details{
  This is a method for the generic function \code{\link[stats]{update}}
  for the class \code{"symbolmap"} of graphics symbol maps.
  It updates the \code{object} using the parameters given in the
  extra arguments \code{\dots}.

  The extra arguments must be given in
  the form \code{name=value} and must be recognisable to
  \code{\link{symbolmap}}. They override any parameters
  of the same name in \code{object}.
}
\value{
  Another object of class \code{"symbolmap"}.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/},
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak.
}
\seealso{
  \code{\link{symbolmap}} to create a graphics symbol map.
}
\examples{
  g <- symbolmap(size=function(x) x/50)
  g
  update(g, range=c(0,1))
  update(g, size=42)
  update(g, shape="squares", range=c(0,1))
}
\keyword{spatial}
\keyword{hplot}
