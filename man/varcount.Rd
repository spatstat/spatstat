\name{varcount}
\alias{varcount}
\title{
  Predicted Variance of the Number of Points 
}
\description{
  Given a fitted point process model, calculate the predicted variance
  of the number of points in a nominated set \code{B}.
}
\usage{
varcount(model, B, \dots, dimyx = NULL)
}
\arguments{
  \item{model}{
    A fitted point process model
    (object of class \code{"ppm"}, \code{"kppm"} or \code{"dppm"}).
  }
  \item{B}{
    A window (object of class \code{"owin"} specifying the region in
    which the points are counted.
    Alternatively a pixel image (object of class \code{"im"})
    or a function of spatial coordinates specifying a numerical weight
    for each random point.
  }
  \item{\dots}{
    Additional arguments passed to \code{B} when it is a function.
  }
  \item{dimyx}{
    Spatial resolution for the calculations.
    Argument passed to \code{\link{as.mask}}.
  }
}
\details{
  This command calculates the variance of the number of points
  falling in a specified window \code{B} according to the \code{model}.
  It can also calculate the variance of a sum of weights attached
  to each random point.

  The \code{model} should be a fitted point process model
  (object of class \code{"ppm"}, \code{"kppm"} or \code{"dppm"}).

  \itemize{
    \item{
      If \code{B} is a window, this command calculates the variance
      of the number of points falling in \code{B}, according to the
      fitted \code{model}.

      If the \code{model} depends on spatial covariates other than the
      Cartesian coordinates, then \code{B} should be a subset of the
      domain in which these covariates are defined. 
    }
    \item{
      If \code{B} is a pixel image,
      this command calculates the variance of
      \eqn{T = \sum_i B(x_i)}{T = sum[i] B(x[i])},
      the sum of the values of \code{B} over all random points
      falling in the domain of the image.

      If the \code{model} depends on spatial covariates other than the
      Cartesian coordinates, then the domain of the pixel image,
      \code{as.owin(B)}, should be a subset of the domain in which these
      covariates are defined. 
    }
    \item{
      If \code{B} is a \code{function(x,y)} or \code{function(x,y,...)}
      this command calculates the variance of
      \eqn{T = \sum_i B(x_i)}{T = sum[i] B(x[i])},
      the sum of the values of \code{B} over all random points
      falling inside the window \code{W=as.owin(model)}, the window
      in which the original data were observed.
    }
  }

  The variance calculation involves the intensity and the
  pair correlation function of the model.
  The calculation is exact (up to discretisation error)
  for models of class \code{"kppm"} and \code{"dppm"},
  and for Poisson point process models of class \code{"ppm"}.
  For Gibbs point process models of class \code{"ppm"} the
  calculation depends on the Poisson-saddlepoint approximations
  to the intensity and pair correlation function, which are rough
  approximations. The approximation is not yet implemented
  for some Gibbs models.
}
\value{
  A single number.
}
\author{
  \spatstatAuthors
}
\seealso{
  \code{\link{predict.ppm}},
  \code{\link{predict.kppm}},
  \code{\link{predict.dppm}}
}
\examples{
   fitT <- kppm(redwood ~ 1, "Thomas")
   B <- owin(c(0, 0.5), c(-0.5, 0))
   varcount(fitT, B)

   fitS <- ppm(swedishpines ~ 1, Strauss(9))
   BS <- square(50)
   varcount(fitS, BS)
}
\keyword{spatial}
\keyword{models}

