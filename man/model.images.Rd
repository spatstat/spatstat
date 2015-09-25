\name{model.images}  %DontDeclareMethods
\alias{model.images}
\alias{model.images.ppm}
\alias{model.images.dppm}
\alias{model.images.kppm}
\alias{model.images.lppm}
\alias{model.images.slrm}
\title{Compute Images of Constructed Covariates}
\description{
  For a point process model fitted to spatial point pattern data,
  this function computes pixel images of the covariates
  in the design matrix.
}
\usage{
  model.images(object, ...)

  \method{model.images}{ppm}(object, W = as.owin(object), ...)

  \method{model.images}{kppm}(object, W = as.owin(object), ...)

  \method{model.images}{dppm}(object, W = as.owin(object), ...)

  \method{model.images}{lppm}(object, L = as.linnet(object), ...)

  \method{model.images}{slrm}(object, ...)
}
\arguments{
  \item{object}{
    The fitted point process model. An object of class \code{"ppm"}
    or \code{"kppm"} or \code{"lppm"} or \code{"slrm"} or \code{"dppm"}.
  }
  \item{W}{
    A window (object of class \code{"owin"}) in which the
    images should be computed. Defaults to the window
    in which the model was fitted.
  }
  \item{L}{
    A linear network (object of class \code{"linnet"}) in which the
    images should be computed. Defaults to the network
    in which the model was fitted.
  }
  \item{\dots}{
    Other arguments (such as \code{na.action}) passed to
    \code{\link{model.matrix.lm}}.
  }
}
\details{
  This command is similar to \code{\link{model.matrix.ppm}} except
  that it computes pixel images of the covariates,
  instead of computing the covariate values at certain points only.

  The \code{object} must be a fitted spatial point process model
  object of class \code{"ppm"} (produced by the model-fitting
  function \code{\link{ppm}}) or class \code{"kppm"} (produced by the
  fitting function \code{\link{kppm}})
  or class \code{"dppm"} (produced by the
  fitting function \code{\link{dppm}}) or class \code{"lppm"} (produced
  by \code{\link{lppm}}) or class \code{"slrm"} (produced by
  \code{\link{slrm}}). 

  The spatial covariates required by the model-fitting procedure
  are computed at every pixel location in the window \code{W}.
  For \code{lppm} objects, the covariates are computed at every
  location on the network \code{L}. For \code{slrm} objects, the
  covariates are computed on the pixels that were used to fit the
  model.

  Note that the spatial covariates computed here
  are not the original covariates that were supplied when fitting the
  model. Rather, they are the covariates that actually appear in the
  loglinear representation of the (conditional) intensity
  and in the columns of the design matrix. For example, they might include
  dummy or indicator variables for different levels of a factor,
  depending on the contrasts that are in force.

  The pixel resolution is determined by \code{W} 
  if \code{W} is a mask (that is \code{W$type = "mask"}).
  Otherwise, the pixel resolution is determined by
  \code{\link{spatstat.options}}.

  The format of the result depends on whether the original point pattern
  data were marked or unmarked.
  \itemize{
    \item
    If the original dataset was unmarked,
    the result is a named list of pixel images (objects of class
    \code{"im"}) containing the values of the spatial covariates.
    The names of the list elements are the names of the covariates
    determined by \code{\link{model.matrix.lm}}.
    The result is also of class \code{"solist"} so that it can
    be plotted immediately.
    \item 
    If the original dataset was a multitype point pattern,
    the result is a \code{\link{hyperframe}}
    with one column for each possible type of points.
    Each column is a named list of pixel images (objects of class
    \code{"im"}) containing the values of the spatial covariates.
    The row names of the hyperframe are the names of the covariates
    determined by \code{\link{model.matrix.lm}}.
  }
}
\value{
  A list (of class \code{"solist"}) or
  array (of class \code{"hyperframe"}) containing
  pixel images (objects of class \code{"im"}).
  For \code{model.images.lppm}, the images are also of class \code{"linim"}.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\seealso{
  \code{\link{model.matrix.ppm}},
  \code{\link[stats]{model.matrix}},
  \code{\link{ppm}},
  \code{\link{ppm.object}},
  \code{\link{lppm}},
  \code{\link{dppm}},
  \code{\link{kppm}},
  \code{\link{slrm}},
  \code{\link{im}},
  \code{\link{im.object}},
  \code{\link{plot.solist}},
  \code{\link{spatstat.options}}
}

\examples{
   fit <- ppm(cells ~ x)
   model.images(fit)
   B <- owin(c(0.2, 0.4), c(0.3, 0.8))
   model.images(fit, B)
   fit2 <- ppm(cells ~ cut(x,3))
   model.images(fit2)
   fit3 <- slrm(japanesepines ~ x)
   model.images(fit3)
   fit4 <- ppm(amacrine ~ marks + x)
   model.images(fit4)
}
\keyword{spatial}
\keyword{models}
