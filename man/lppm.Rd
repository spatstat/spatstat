\name{lppm}
\alias{lppm}
\title{
  Fit Point Process Model to Point Pattern on Linear Network
}
\description{
  Fit a point process model to a point pattern dataset on a linear network
}
\usage{
lppm(X, ..., eps=NULL, nd=1000)
}
\arguments{
  \item{X}{
    Object of class \code{"lpp"} specifying a point pattern
    on a linear network.
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat]{ppm}}.
  }
  \item{eps}{
    Optional. Spacing between dummy points along each segment of the
    network. 
  }
  \item{nd}{
    Optional. Number of dummy points equally spaced along each segment
    of the network. Ignored if \code{eps} is given.
  }
}
\details{
  This function fits a point process model to data that specify
  a point pattern on a linear network. It is a counterpart of
  the model-fitting function \code{\link[spatstat]{ppm}} designed
  to work with objects of class \code{"lpp"} instead of \code{"ppp"}.

  The argument \code{X} should be an object of class \code{"lpp"}
  (created by the command \code{\link{lpp}}) specifying a point pattern
  on a linear network.

  The subsequent arguments \code{...} will be passed to
  \code{\link[spatstat]{ppm}}. They specify the form of the model.
}
\value{
  An object of class \code{"lppm"} representing the fitted model.
  There are methods for \code{print}, \code{predict},
  \code{coef} and similar functions.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\seealso{
  \code{\link{methods.lppm}},
  \code{\link{predict.lppm}},
  \code{\link{ppm}},
  \code{\link{lpp}}.
}
\examples{
  example(lpp)
  lppm(X, ~1)
  lppm(X, ~x)
}
\references{
  Ang, Q.W. (2010)
  \emph{Statistical methodology for events on a network}.
  Master's thesis, School of Mathematics and Statistics, University of
  Western Australia.
  
  Ang, Q.W., Baddeley, A. and Nair, G. (2012)
  Geometrically corrected second-order analysis of 
  events on a linear network, with applications to
  ecology and criminology.
  To appear in \emph{Scandinavian Journal of Statistics}.

  McSwiggan, G., Nair, M.G. and Baddeley, A. (2012)
  Fitting Poisson point process models to events 
  on a linear network. Manuscript in preparation.
}
\keyword{spatial}
\keyword{models}
