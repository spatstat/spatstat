\name{plot.mppm}
\alias{plot.mppm}
\title{plot a Fitted Multiple Point Process Model}
\description{
  Given a point process model fitted to multiple point patterns
  by \code{\link{mppm}},
  compute spatial trend or conditional intensity surface of the model,
  in a form suitable for plotting, and (optionally) plot this
  surface.
}
\usage{
  \method{plot}{mppm}(x, ..., trend = TRUE, cif = FALSE,
                  how="image")
}
\arguments{
  \item{x}{
    A point process model fitted to multiple point patterns,
    typically obtained from
    the model-fitting algorithm \code{\link{mppm}}.
    An object of class \code{"mppm"}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{plot.ppm}} or
    \code{\link{plot.listof}} controlling the plot.
  }
  \item{trend}{
    logical flag; if \code{TRUE}, the spatial trend surface will be produced.
  }
  \item{cif}{
    logical flag; if \code{TRUE}, the conditional intensity surface will be
    produced.
  }
  \item{how}{
    character string vector indicating the style of
    plot to be performed. 
  }
}
\value{
  \code{NULL}.
}
\details{
  This is the \code{plot} method for the class \code{"mppm"}
  of point process models fitted to multiple point patterns
  (see \code{\link{mppm}}).

  It invokes \code{\link{subfits}} to compute the fitted model for
  each individual point pattern dataset, then calls
  \code{\link{plot.ppm}} to plot these individual models. These
  individual plots are displayed using \code{\link{plot.listof}},
  which generates either a series of separate plot frames or an
  array of plot panels on a single page. 
}
\seealso{
  \code{\link{plot.ppm}},
  \code{\link{mppm}},
  \code{\link{plot.listof}}
}
\references{
  Baddeley, A., Rubak, E. and Turner, R. (2015)
  \emph{Spatial Point Patterns: Methodology and Applications with R}.
  London: Chapman and Hall/CRC Press. 
}
\author{
  Adrian Baddeley, Ida-Maria Sintorn and Leanne Bischoff.
  Implemented by 
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{hplot}
\keyword{models}
