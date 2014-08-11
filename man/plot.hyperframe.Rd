\name{plot.hyperframe}
\alias{plot.hyperframe}
\title{Plot Entries in a Hyperframe}
\description{
  Plots the entries in a hyperframe, in a series of
  panels, one panel for each row of the hyperframe.
}
\usage{
   \method{plot}{hyperframe}(x, e, ..., main, arrange=TRUE,
                            nrows=NULL, ncols=NULL,
                            parargs=list(mar=c(1,1,3,1) * marsize),
                            marsize=0.1)
}
\arguments{
  \item{x}{
    Data to be plotted. A hyperframe (object of class \code{"hyperframe"}, see
    \code{\link{hyperframe}}).
  }
  \item{e}{
    How to plot each row. Optional. \R language expression that will be
    evaluated in each row of the hyperframe to generate the plots.
  }
  \item{\dots}{
    Extra arguments controlling the plot (when \code{e} is missing).
  }
  \item{main}{Overall title for the array of plots.}
  \item{arrange}{
    Logical flag indicating whether to plot the objects
    side-by-side on a single page (\code{arrange=TRUE})
    or plot them individually in a succession of frames
    (\code{arrange=FALSE}).
  }
  \item{nrows,ncols}{
    Optional. The number of rows/columns in the plot layout
    (assuming \code{arrange=TRUE}).
    You can specify either or both of these numbers.
  }
  \item{parargs}{
    Optional list of arguments passed to \code{\link{par}} before
    plotting each panel. Can be used to control margin sizes, etc.
  }
  \item{marsize}{
    Optional scale parameter controlling the sizes of margins between
    the panels.
  }
}
\details{
  This is the \code{plot} method for the class \code{"hyperframe"}.

  The argument \code{x} must be a hyperframe (like a data frame,
  except that the entries can be objects of any class; see
  \code{\link{hyperframe}}).

  This function generates a series of plots, one plot for each
  row of the hyperframe. If \code{arrange=TRUE} (the default), then
  these plots are arranged in a neat array of panels within a single
  plot frame. If \code{arrange=FALSE}, the plots are simply executed
  one after another.

  Exactly what is plotted, and how it is plotted, depends on the
  argument \code{e}. The default (if \code{e} is missing) is to plot
  only the first column of \code{x}. Each entry in the first column
  is plotted using the generic \code{\link{plot}} command, together with
  any extra arguments given in \code{\dots}.

  If \code{e} is present, it should be an \R language expression
  involving the column names of \code{x}.
  The expression will be evaluated once for each row of \code{x}.
  It will be evaluated in an environment where each column name of
  \code{x} is interpreted as meaning the object in that column
  in the current row.
  See the Examples.
}
\value{
  \code{NULL}.
}
\seealso{
  \code{\link{hyperframe}},
  \code{\link{with.hyperframe}}
}
\examples{
   H <- hyperframe(id=1:10)
   H$X <- with(H, rpoispp(100))
   H$D <- with(H, distmap(X))
   # points only
   plot(H[,"X"])
   plot(H, plot(X, main=id))
   # points superimposed on images
   plot(H, {plot(D, main=id); plot(X, add=TRUE)})
   
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{hplot}
