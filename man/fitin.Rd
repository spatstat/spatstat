\name{fitin.ppm}
\alias{fitin}
\alias{fitin.ppm}
\alias{fitin.profilepl}
\title{Extract the Interaction from a Fitted Point Process Model}
\description{
  Given a point process model that has been fitted to point pattern
  data, this function extracts the interpoint interaction part of the
  model as a separate object.
}
\usage{
fitin(object)

\method{fitin}{ppm}(object)

\method{fitin}{profilepl}(object)
}
\arguments{
  \item{object}{A fitted point process model (object of class
    \code{"ppm"} or \code{"profilepl"}).
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

  For a pairwise interaction, the plot of the fitted interaction
  shows the pair interaction function (the contribution to the
  probability density from a pair of points as a function of the
  distance between them). For a higher-order interaction, the plot shows
  the strongest interaction (the value most different from 1)
  that could ever arise at the given distance. 

  The fitted interaction coefficients can also be extracted
  from this object using \code{\link{coef}}.
}
\value{
  An object of class \code{"fii"} representing the fitted
  interpoint interaction. This object can be printed and plotted.
}
\author{
  \spatstatAuthors.
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
   model <- ppm(swedishpines ~1, PairPiece(seq(3,19,by=4)))
   f <- fitin(model)
   f
   plot(f)

# extract fitted interaction coefficients
   coef(f)

   # multitype
   # fit the stationary multitype Strauss process to `amacrine'
   r <- 0.02 * matrix(c(1,2,2,1), nrow=2,ncol=2)
   model <- ppm(amacrine ~1, MultiStrauss(r))
   f <- fitin(model)
   f
   plot(f)
}
\keyword{spatial}
\keyword{models}
