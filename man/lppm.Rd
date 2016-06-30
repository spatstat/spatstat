\name{lppm}
\alias{lppm}
\alias{lppm.formula}
\alias{lppm.lpp}
\title{
  Fit Point Process Model to Point Pattern on Linear Network
}
\description{
  Fit a point process model to a point pattern dataset on a linear network
}
\usage{
lppm(X, ...)

\method{lppm}{formula}(X, interaction=NULL, ..., data=NULL)

\method{lppm}{lpp}(X, ..., eps=NULL, nd=1000)
}
\arguments{
  \item{X}{
    Either an object of class \code{"lpp"} specifying a point pattern
    on a linear network, or a \code{formula} specifying the
    point process model.
  }
  \item{\dots}{
    Arguments passed to \code{\link{ppm}}.
  }
  \item{interaction}{
    An object of class \code{"interact"}
    describing the point process interaction
    structure, or \code{NULL} indicating that a Poisson process (stationary
    or nonstationary) should be fitted.
  }
  \item{data}{
    Optional. The values of spatial covariates (other than the Cartesian
    coordinates) required by the model.
    A list whose entries are images,
    functions, windows, tessellations or single numbers. 
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
  the model-fitting function \code{\link{ppm}} designed
  to work with objects of class \code{"lpp"} instead of \code{"ppp"}.

  The function \code{lppm} is generic, with methods for
  the classes \code{formula} and \code{lppp}.

  In \code{lppm.lpp}
  the first argument \code{X} should be an object of class \code{"lpp"}
  (created by the command \code{\link{lpp}}) specifying a point pattern
  on a linear network.

  In \code{lppm.formula},
  the first argument is a \code{formula} in the \R language
  describing the spatial trend model to be fitted. It has the general form
  \code{pattern ~ trend} where the left hand side \code{pattern} is usually
  the name of a point pattern on a linear network
  (object of class \code{"lpp"})
  to which the model should be fitted, or an expression which evaluates
  to such a point pattern;
  and the right hand side \code{trend} is an expression specifying the
  spatial trend of the model.

  Other arguments \code{...} are passed from \code{lppm.formula}
  to \code{lppm.lpp} and from \code{lppm.lpp} to \code{\link{ppm}}.
}
\value{
  An object of class \code{"lppm"} representing the fitted model.
  There are methods for \code{print}, \code{predict},
  \code{coef} and similar functions.
}
\author{
  \adrian
  
  
}
\seealso{
  \code{\link{methods.lppm}},
  \code{\link{predict.lppm}},
  \code{\link{ppm}},
  \code{\link{lpp}}.
}
\examples{
  X <- runiflpp(15, simplenet)
  lppm(X ~1)
  lppm(X ~x)
  marks(X) <- factor(rep(letters[1:3], 5))
  lppm(X ~ marks)
  lppm(X ~ marks * x)
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
  \emph{Scandinavian Journal of Statistics} \bold{39}, 591--617.

  McSwiggan, G., Nair, M.G. and Baddeley, A. (2012)
  Fitting Poisson point process models to events 
  on a linear network. Manuscript in preparation.
}
\keyword{spatial}
\keyword{models}
