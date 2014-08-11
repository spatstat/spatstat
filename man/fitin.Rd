\name{fitin.ppm}
\alias{fitin}
\alias{fitin.ppm}
\title{Extract the Interaction from a Fitted Point Process Model}
\description{
  Given a point process model that has been fitted to point pattern
  data, this function extracts the interpoint interaction part of the
  model as a separate object.
}
\usage{
fitin(object)
\method{fitin}{ppm}(object)
}
\arguments{
  \item{object}{A fitted point process model (object of class
    \code{"ppm"}).
  }
}
\details{
  An object of class \code{"ppm"} describes a fitted point process
  model. It contains information about the original data to which the
  model was fitted, the spatial trend that was fitted, the
  interpoint interaction that was fitted, and other data.
  See \code{\link{ppm.object}}) for details of this class.

  The function \code{fitin} extracts from this model the information about the
  fitted interpoint interaction only.
  The information is organised as an object of class \code{"fii"}
  (fitted interpoint interaction).
  This object can be printed or plotted.

  Users may find this a convenient way to plot the
  fitted interpoint interaction term, as shown in the Examples.

  The fitted interaction coefficients can also be extracted
  from this object using \code{\link{coef}}.
}
\value{
  An object of class \code{"fii"} representing the fitted
  interpoint interaction. This object can be printed and plotted.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  Methods for handling fitted interactions: 
  \code{\link{methods.fii}}, \code{\link{reach.fii}},
  \code{\link{as.interact.fii}}.

  Background:
  \code{\link{ppm}},
  \code{\link{ppm.object}}.
}
\examples{
   # unmarked 
   model <- ppm(swedishpines, ~1, PairPiece(seq(3,19,by=4)))
   f <- fitin(model)
   f
   plot(f)

# extract fitted interaction coefficients
   coef(f)

   # multitype
   # fit the stationary multitype Strauss process to `amacrine'
   r <- 0.02 * matrix(c(1,2,2,1), nrow=2,ncol=2)
   model <- ppm(amacrine, ~1, MultiStrauss(c("off","on"), r))
   f <- fitin(model)
   f
   plot(f)
}
\keyword{spatial}
\keyword{models}
