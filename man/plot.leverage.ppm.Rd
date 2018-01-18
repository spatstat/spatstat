\name{plot.leverage.ppm}
\alias{plot.leverage.ppm}
\alias{contour.leverage.ppm}
\alias{persp.leverage.ppm}
\title{
  Plot Leverage Function
}
\description{
  Generate a pixel image plot, or a contour plot, or a perspective plot,
  of a leverage function that has been computed by \code{\link{leverage.ppm}}.
}
\usage{
 \method{plot}{leverage.ppm}(x, \dots,
                             what=c("smooth", "nearest", "exact"),
                             showcut=TRUE,
                             args.cut=list(drawlabels=FALSE), 
                             multiplot=TRUE)

 \method{contour}{leverage.ppm}(x, \dots,
                             what=c("smooth", "nearest"),
                             showcut=TRUE,
                             args.cut=list(col=3, lwd=3, drawlabels=FALSE), 
                             multiplot=TRUE)

 \method{persp}{leverage.ppm}(x, \dots,
               what=c("smooth", "nearest"),
               main, zlab="leverage")
}
\arguments{
  \item{x}{
    Leverage function (object of class \code{"leverage.ppm"}) computed by
    \code{\link{leverage.ppm}}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{plot.im}} or \code{\link{contour.im}}
    or \code{\link{persp.im}} controlling the plot.
  }
  \item{what}{
    Character string (partially matched)
    specifying the values to be plotted. See Details.
  }
  \item{showcut}{
    Logical. If \code{TRUE}, a contour line is plotted at the
    level equal to the theoretical mean of the leverage.
  }
  \item{args.cut}{
    Optional list of arguments passed to
    \code{\link[graphics]{contour.default}} to control the
    plotting of the contour line for the mean leverage.
  }
  \item{multiplot}{
    Logical value indicating whether it is permissible to display
    several plot panels.
  }
  \item{main}{
    Optional main title. A character string or character vector.
  }
  \item{zlab}{
    Label for the \eqn{z} axis. A character string.
  }
}
\details{
  These functions are the \code{plot}, \code{contour} and \code{persp} methods
  for objects of class \code{"leverage.ppm"}.
  Such objects are computed by the command \code{\link{leverage.ppm}}.

  The \code{plot} method displays the leverage function 
  as a colour pixel image using \code{\link{plot.im}},
  and draws a single contour line at the mean leverage value
  using \code{\link{contour.default}}.
  Use the argument \code{clipwin} to restrict the plot to a subset
  of the full data.

  The \code{contour} method displays the leverage function as a contour
  plot, and also draws a single contour line at the mean leverage value,
  using \code{\link{contour.im}}.
  
  The \code{persp} method displays the leverage function as a surface
  in perspective view, using \code{\link{persp.im}}.

  Since the exact values of leverage are computed only
  at a finite set of quadrature locations, there are several options for
  these plots:
  \describe{
    \item{\code{what="smooth"}:}{
      (the default) an image plot showing a smooth function,
      obtained by applying kernel smoothing to the exact leverage values;
    }
    \item{\code{what="nearest"}:}{
      an image plot showing a piecewise-constant function,
      obtained by taking the exact leverage value at the nearest
      quadrature point;
    }
    \item{\code{what="exact"}:}{
      a symbol plot showing the exact values of leverage as
      circles, centred at the quadrature points,
      with diameters proportional to leverage.
    }
  }
  The pixel images are already contained in the object \code{x}
  and were computed by \code{\link{leverage.ppm}};
  the resolution of these images is controlled
  by arguments to \code{\link{leverage.ppm}}.
}
\value{
  Same as for \code{\link{plot.im}}, \code{\link{contour.im}}
  and \code{\link{persp.im}} respectively.
}
\references{
  Baddeley, A., Chang, Y.M. and Song, Y. (2013)
  Leverage and influence diagnostics for spatial point process models.
  \emph{Scandinavian Journal of Statistics} \bold{40}, 86--104.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{leverage.ppm}}.
}
\examples{
   X <- rpoispp(function(x,y) { exp(3+3*x) })
   fit <- ppm(X ~x+y)
   lef <- leverage(fit)
   plot(lef)
   contour(lef)
   persp(lef)
}
\keyword{spatial}
\keyword{models}

