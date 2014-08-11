\name{as.interact}
\alias{as.interact}
\alias{as.interact.fii}
\alias{as.interact.interact}
\alias{as.interact.ppm}
\title{Extract Interaction Structure}
\description{
  Extracts the interpoint interaction structure from
  a point pattern model.
}
\usage{
as.interact(object)
\method{as.interact}{fii}(object)
\method{as.interact}{interact}(object)
\method{as.interact}{ppm}(object)
}
\arguments{
  \item{object}{A fitted point process model (object of class
    \code{"ppm"}) or an interpoint interaction structure
    (object of class \code{"interact"}).
  }
}
\details{
  The function \code{as.interact} extracts 
  the interpoint interaction structure from a suitable object.
  
  An object of class \code{"interact"} describes an interpoint
  interaction structure, before it has been fitted to point pattern
  data. The irregular parameters of the interaction (such as the
  interaction range) are fixed, but the regular parameters
  (such as interaction strength) are undetermined.
  Objects of this class are created by the functions
  \code{\link{Poisson}}, \code{\link{Strauss}} and so on.
  The main use of such objects is in a call to \code{\link{ppm}}.
  
  The function \code{as.interact} is generic, with methods for the classes
  \code{"ppm"}, \code{"fii"} and \code{"interact"}. 
  The result is an object of class \code{"interact"} which can be printed.
}
\section{Note on parameters}{
  This function does \bold{not} extract the fitted coefficients
  of the interaction. To extract the fitted interaction including
  the fitted coefficients, use \code{\link{fitin}}.
}
\value{
  An object of class \code{"interact"} representing the 
  interpoint interaction. This object can be printed and plotted.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{fitin}},
  \code{\link{ppm}}.
}
\examples{
   data(cells)
   model <- ppm(cells, ~1, Strauss(0.07))
   f <- as.interact(model)
   f
}
\keyword{spatial}
\keyword{models}
