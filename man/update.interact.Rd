\name{update.interact}
\alias{update.interact}
\title{
  Update an Interpoint Interaction
}
\description{
  This command updates the \code{object} using the
  arguments given.
}
\usage{
\method{update}{interact}(object, \dots)
}
\arguments{
  \item{object}{
    Interpoint interaction (object of class \code{"interact"}).
  }
  \item{\dots}{
    Additional or replacement values of parameters of \code{object}.
  }
}
\details{
  This is a method for the generic function \code{\link[stats]{update}}
  for the class \code{"interact"} of interpoint interactions.
  It updates the \code{object} using the parameters given in the
  extra arguments \code{\dots}.

  The extra arguments must be given in
  the form \code{name=value} and must be recognisable to
  the interaction object. They override any parameters
  of the same name in \code{object}.
}
\value{
  Another object of class \code{"interact"}, equivalent to \code{object}
  except for changes in parameter values.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au},
 Rolf Turner \email{r.turner@auckland.ac.nz}
 and Ege Rubak \email{rubak@math.aau.dk}.
}
\seealso{
  \code{\link{update.ppm}}
}
\examples{
  Str <- Strauss(r=1)
  Str
  update(Str, r=2)

  M <- MultiStrauss(radii=matrix(1,2,2))
  update(M, types=c("on", "off"))
}
\keyword{spatial}
\keyword{models}
