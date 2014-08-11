\name{as.fv} %DontDeclareMethods
\alias{as.fv}
\alias{as.fv.fv}
\alias{as.fv.fasp}
\alias{as.fv.data.frame}
\alias{as.fv.matrix}
\alias{as.fv.minconfit}
\alias{as.fv.kppm}
\alias{as.fv.bw.optim}
\title{Convert Data To Class fv}
\description{
  Converts data into a function table (an object of class \code{"fv"}).
}
\usage{
  as.fv(x)

  \method{as.fv}{fv}(x)

  \method{as.fv}{data.frame}(x)

  \method{as.fv}{matrix}(x)

  \method{as.fv}{fasp}(x)

  \method{as.fv}{minconfit}(x)

  \method{as.fv}{kppm}(x)

  \method{as.fv}{bw.optim}(x)
}
\arguments{
  \item{x}{Data which will be converted into a function table}
}
\value{
  An object of class \code{"fv"} (see \code{\link{fv.object}}).
}
\details{
  This command converts data \code{x}, that
  could be interpreted as the values of a function,
  into a function value table (object of the class \code{"fv"}
  as described in \code{\link{fv.object}}). This object can then
  be plotted easily using \code{\link{plot.fv}}.

  The dataset \code{x} may be any of the following:
  \itemize{
    \item
    an object of class \code{"fv"};
    \item
    a matrix or data frame with at least two columns;
    \item
    an object of class \code{"fasp"}, representing an array of
    \code{"fv"} objects.
    \item
    an object of class \code{"minconfit"}, giving the results
    of a minimum contrast fit by the command \code{\link{mincontrast}}.
    The 
    \item
    an object of class \code{"kppm"}, representing a fitted
    Cox or cluster point process model, obtained from the
    model-fitting command \code{\link{kppm}};
    \item
    an object of class \code{"bw.optim"}, representing an optimal
    choice of smoothing bandwidth by a cross-validation method, obtained
    from commands like \code{\link{bw.diggle}}.
  }
  The function \code{as.fv} is generic, with methods for each of the
  classes listed above. The behaviour is as follows:
  \itemize{
    \item
    If \code{x} is an object of class \code{"fv"}, it is
    returned unchanged.
    \item
    If \code{x} is a matrix or data frame, 
    the first column is interpreted
    as the function argument, and subsequent columns are interpreted as
    values of the function computed by different methods.
    \item
    If \code{x} is an object of class \code{"fasp"}
    representing an array of \code{"fv"} objects, 
    these are combined into a single \code{"fv"} object.
    \item
    If \code{x} is an object of class \code{"minconfit"},
    or an object of class \code{"kppm"} that was fitted by minimum contrast,
    the result is a function table containing the
    observed summary function and the best fit summary function.
    \item
    If \code{x} is an object of class \code{"bw.optim"},
    the result is a function table of the optimisation criterion
    as a function of the smoothing bandwidth.
  }
}
\examples{
  r <- seq(0, 1, length=101)
  x <- data.frame(r=r, y=r^2)
  as.fv(x)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
