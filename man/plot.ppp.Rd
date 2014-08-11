\name{plot.ppp}
\alias{plot.ppp}
\title{plot a Spatial Point Pattern}
\description{
  Plot a two-dimensional spatial point pattern
}
\usage{
 \method{plot}{ppp}(x, main, \dots, chars=NULL, cols=NULL,
use.marks=TRUE, which.marks=1, add=FALSE,
maxsize=NULL, markscale=NULL,zap=0.01)
}
\arguments{
  \item{x}{
    The spatial point pattern to be plotted.
    An object of class \code{"ppp"},
    or data which can be converted into 
    this format by \code{\link{as.ppp}()}.
  }
  \item{main}{
    text to be displayed as a title above the plot.
  }
  \item{\dots}{
    extra arguments that will be passed to the plotting functions
    \code{\link{plot.default}}, \code{\link{points}} and/or
    \code{\link{symbols}}.
  }
  \item{chars}{
    plotting character(s) used to plot points.
  }
  \item{cols}{
    the colour(s) used to plot points.
  }
  \item{use.marks}{
    logical flag; if \code{TRUE}, plot points using a different
    plotting symbol for each mark;
    if \code{FALSE}, only the locations of the points will be plotted,
    using \code{\link{points}()}.
  }
  \item{which.marks}{
    Index determining which column of marks to use,
    if the marks of \code{x} are a data frame.
    A character string or an integer.
    Defaults to \code{1} indicating the first column of marks.
  }
  \item{add}{
    logical flag; if \code{TRUE},
    just the points are plotted, over the existing plot.
    A new plot is not created, and
    the window is not plotted.
  }
  \item{maxsize}{
    maximum size of the circles/squares plotted
    when \code{x} is a marked point pattern with 
    numerical marks.
    Incompatible with \code{markscale}.
  }
  \item{markscale}{
    physical scale factor determining the sizes of the
    circles/squares plotted when \code{x} is a marked point pattern with 
    numerical marks. Incompatible with \code{maxsize}.
  }
  \item{zap}{
    Fraction between 0 and 1.
    When \code{x} is a marked point pattern with numerical marks,
    \code{zap} is the smallest mark value
    (expressed as a fraction of the maximum possible mark) that will
    be plotted. 
    Any points which have marks smaller in absolute value
    than \code{zap * max(abs(marks(x)))} will not be plotted.
  }
}
\value{
  \code{NULL}, or a vector giving the correspondence between 
  mark values and plotting characters.
}
\details{
  This is the \code{plot} method for 
  point pattern datasets (of class \code{"ppp"}, see \code{\link{ppp.object}}).

  First the observation window \code{x$window} is plotted.
  Then the points themselves are plotted,
  in a fashion that depends on their marks,
  as follows.
  \describe{
    \item{unmarked point pattern:}{
      If the point pattern does not have marks, or if \code{use.marks = FALSE},
      then the locations of all points will be plotted 
      using a single plot character
    }
    \item{multitype point pattern:}{
      If \code{x$marks} is a factor, then 
      each level of the factor is 
      represented by a different plot character.
    }
    \item{continuous marks:}{
      If \code{x$marks} is a numeric vector,
      the marks are rescaled to the unit interval and
      each point is represented by a circle 
      with radius proportional to the rescaled mark
      (if the value is positive) or a square with side proportional
      to the absolute value of the rescaled mark (if the value is negative).
    }
    \item{other kinds of marks:}{
      If \code{x$marks} is neither numeric nor a factor,
      then each possible mark will be represented by a
      different plotting character. The default is to
      represent the \eqn{i}th smallest mark value by
      \code{points(..., pch=i)}.
    }
  }

  Plotting of the window \code{x$window} is performed by
  \code{\link{plot.owin}}. This plot may be modified 
  through the \code{...} arguments. In particular the
  extra argument \code{border} determines
  the colour of the window.

  Plotting of the points themselves is performed
  by the function \code{\link{points}}, except for the case of
  continuous marks, where it is performed by \code{\link{symbols}}.
  Their plotting behaviour may be modified through the \code{...}
  arguments.

  The argument \code{chars} determines the plotting character
  or characters used to display the points (in all cases except
  for the case of continuous marks). For an unmarked point pattern,
  this should be a single integer or character determining a
  plotting character (see \code{par("pch")}).
  For a multitype point pattern, \code{chars} should be a vector
  of integers or characters, of the same length
  as \code{levels(x$marks)}, and then the \eqn{i}th level or type 
  will be plotted using character \code{chars[i]}.

  If \code{chars} is absent, but there is an extra argument
  \code{pch}, then this will determine the plotting character for
  all points.
  
  The argument \code{cols} determines the colour or colours used to
  display the points. For an unmarked point pattern, or a
  marked point pattern with continuous marks, this should be a character string
  determining a colour. For a multitype point pattern, \code{cols}
  should be a character vector, of the same length
  as \code{levels(x$marks)}. The \eqn{i}th level or type will
  be plotted using colour \code{cols[i]}.

  If \code{cols} is absent, the colour used to plot \emph{all} the
  points may be determined by the extra argument \code{fg}
  (for multitype point patterns) or the extra argument \code{col}
  (for all other cases). Note that \code{col} will also reset the
  colour of the window.

  The arguments \code{maxsize} and \code{markscale} 
  incompatible. They control the physical size of the circles and
  squares which represent the marks in a point pattern with continuous
  marks. If \code{markscale} is given, then a mark value of \code{m}
  is plotted as a circle of radius \code{m * markscale}
  (if \code{m} is positive) or a square of side \code{abs(m) * markscale}
  (if \code{m} is negative). If \code{maxsize} is given, then the
  largest mark in absolute value, \code{mmax=max(abs(x$marks))},
  will be scaled to have physical size \code{maxsize}.

  The user can set the default values of these plotting parameters
  using \code{\link{spatstat.options}("par.points")}.
  
  To zoom in (to view only a subset of the point pattern at higher
  magnification), use the graphical arguments
  \code{xlim} and \code{ylim} to specify the rectangular field of view.

  The value returned by this plot function can be used to make a
  suitable legend, as shown in the examples.
}
\section{Removing White Space Around The Plot}{
  A frequently-asked question is: How do I remove the white space around
  the plot? Currently \code{plot.ppp} uses the base graphics system of
  \R, so the space around the plot is controlled by parameters 
  to \code{\link{par}}. To reduce the white space, change the
  parameter \code{mar}. Typically, \code{par(mar=rep(0.5, 4))} is
  adequate, if there are no annotations or titles outside the window.
}
\seealso{
  \code{\link{iplot}},
  \code{\link{ppp.object}},
  \code{\link{plot}},
  \code{\link{par}},
  \code{\link{points}},
  \code{\link{plot.owin}},
  \code{\link{symbols}}
}
\examples{
   data(cells)
   plot(cells)

   plot(cells, pch=16)

   # make the plotting symbols larger (for publication at reduced scale)
   plot(cells, cex=2)

   # set it in spatstat.options
   oldopt <- spatstat.options(par.points=list(cex=2))
   plot(cells)
   spatstat.options(oldopt)

   # multitype 
   data(lansing)
   plot(lansing)

   # marked by a real number
   data(longleaf)
   plot(longleaf)

   # just plot the points
   plot(longleaf, use.marks=FALSE)
   plot(unmark(longleaf)) # equivalent

   # controlling COLOURS of points
   plot(cells, cols="blue")
   plot(lansing, cols=c("black", "yellow", "green", 
                        "blue","red","pink"))
   plot(longleaf, fg="blue")

   # make window purple
   plot(lansing, border="purple")
   # make everything purple
   plot(lansing, border="purple", cols="purple", col.main="purple")
 
   # controlling PLOT CHARACTERS
   plot(lansing, chars = 11:16)
   plot(lansing, chars = c("o","h","m",".","o","o"))

   # controlling MARK SCALE
   plot(longleaf, markscale=0.1)

   # draw circles of DIAMETER equal to nearest neighbour distance
   plot(cells \%mark\% nndist(cells), markscale=1/2)

   # making the legend
   data(amacrine)
   v <- plot(amacrine)
   legend(0.2, 1.2, pch=v, legend=names(v))

   # point pattern with multiple marks
   data(finpines)
   plot(finpines, which.marks="height")
}

\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{hplot}
