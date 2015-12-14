\name{plot.ppp}
\alias{plot.ppp}
\title{plot a Spatial Point Pattern}
\description{
  Plot a two-dimensional spatial point pattern
}
\usage{
 \method{plot}{ppp}(x, main, \dots, clipwin=NULL,
                    chars=NULL, cols=NULL,
                    use.marks=TRUE, which.marks=NULL,
                    add=FALSE, type=c("p","n"),
                    legend=TRUE,
                    leg.side=c("left", "bottom", "top", "right"),
                    leg.args=list(),
                    symap=NULL, maxsize=NULL, meansize=NULL, markscale=NULL,
                    zap=0.01,
                    show.window=show.all, show.all=!add, do.plot=TRUE,
                    multiplot=TRUE)
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
  \item{clipwin}{
    Optional. A window (object of class \code{"owin"}).
    Only this subset of the image will be displayed.
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
    A character or integer vector identifying one or more
    columns of marks. 
    If \code{add=FALSE} then
    the default is to plot all columns of marks, in a series of separate
    plots.
    If \code{add=TRUE} then only one column of marks can be plotted,
    and the default is \code{which.marks=1}
    indicating the first column of marks.
  }
  \item{add}{
    logical flag; if \code{TRUE},
    just the points are plotted, over the existing plot.
    A new plot is not created, and
    the window is not plotted.
  }
  \item{type}{
    Type of plot: either \code{"p"} or \code{"n"}.
    If \code{type="p"} (the default), both the points and the observation window
    are plotted. If \code{type="n"}, only the window is plotted.
  }
  \item{legend}{
    Logical value indicating whether to add a legend showing the mapping
    between mark values and graphical symbols (for a marked point pattern).
  }
  \item{leg.side}{
    Position of legend relative to main plot.
  }
  \item{leg.args}{
    List of additional arguments passed to \code{\link{plot.symbolmap}}
    or \code{\link{symbolmap}} to control the legend.
    In addition to arguments documented under
    \code{\link{plot.symbolmap}}, and graphical arguments recognised
    by \code{\link{symbolmap}}, the list may also include the argument
    \code{sep} giving the separation between the main plot and the
    legend, or \code{sep.frac} giving the separation as a fraction
    of the relevant dimension (width or height) of the main plot.
  }
  \item{symap}{
    Optional. The graphical symbol map to be applied to the marks.
    An object of class \code{"symbolmap"}; see
    \code{\link{symbolmap}}.
  }
  \item{maxsize}{
    \emph{Maximum} physical size of the circles/squares plotted
    when \code{x} is a marked point pattern with 
    numerical marks.
    Incompatible with \code{meansize} and \code{markscale}.
    Ignored if \code{symap} is given.
  }
  \item{meansize}{
    \emph{Average} physical size of the circles/squares plotted
    when \code{x} is a marked point pattern with 
    numerical marks.
    Incompatible with \code{maxsize} and \code{markscale}.
    Ignored if \code{symap} is given.
  }
  \item{markscale}{
    physical scale factor determining the sizes of the
    circles/squares plotted when \code{x} is a marked point pattern with 
    numerical marks. Mark value will be multiplied by \code{markscale}
    to determine physical size.
    Incompatible with \code{maxsize} and \code{meansize}.
    Ignored if \code{symap} is given.
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
  \item{show.window}{
    Logical value indicating whether to plot the observation
    window of \code{x}.
  }
  \item{show.all}{
    Logical value indicating whether to plot everything
    including the main title and the observation window of \code{x}.
  }
  \item{do.plot}{
    Logical value determining whether to actually perform the plotting.
  }
  \item{multiplot}{
    Logical value giving permission to display multiple plots.
  }
}
\value{
  (Invisible) object of class \code{"symbolmap"}
  giving the correspondence between 
  mark values and plotting characters.
}
\details{
  This is the \code{plot} method for 
  point pattern datasets (of class \code{"ppp"}, see \code{\link{ppp.object}}).

  First the observation window \code{x$window} is plotted
  (if \code{show.window=TRUE}).
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
      with \emph{diameter} proportional to the rescaled mark
      (if the value is positive) or a square with \emph{side length}
      proportional to the absolute value of the rescaled mark
      (if the value is negative).
    }
    \item{other kinds of marks:}{
      If \code{x$marks} is neither numeric nor a factor,
      then each possible mark will be represented by a
      different plotting character. The default is to
      represent the \eqn{i}th smallest mark value by
      \code{points(..., pch=i)}.
    }
  }

  If there are several columns of marks, and if \code{which.marks} is
  missing or \code{NULL}, then
  \itemize{
    \item 
    if \code{add=FALSE} and \code{multiplot=TRUE}
    the default is to plot all columns of marks, in a series of separate
    plots, placed side-by-side. The plotting is coordinated by
    \code{\link{plot.listof}}, which calls \code{plot.ppp} to make each of
    the individual plots.
    \item
    Otherwise, only one column of marks can be plotted,
    and the default is \code{which.marks=1}
    indicating the first column of marks.
  }

  Plotting of the window \code{Window(x)} is performed by
  \code{\link{plot.owin}}. This plot may be modified 
  through the \code{...} arguments. In particular the
  extra argument \code{border} determines
  the colour of the window, if the window is not a binary mask.

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

  The default colour for the points is a semi-transparent grey,
  if this is supported by the plot device. This behaviour can be
  suppressed (so that the default colour is non-transparent)
  by setting \code{spatstat.options(transparent=FALSE)}.

  The arguments \code{maxsize}, \code{meansize} and \code{markscale} 
  incompatible. They control the physical size of the circles and
  squares which represent the marks in a point pattern with continuous
  marks. The size of a circle is defined as its \emph{diameter};
  the size of a square is its side length.
  If \code{markscale} is given, then a mark value of \code{m}
  is plotted as a circle of diameter \code{m * markscale}
  (if \code{m} is positive) or a square of side \code{abs(m) * markscale}
  (if \code{m} is negative). If \code{maxsize} is given, then the
  largest mark in absolute value, \code{mmax=max(abs(marks(x)))},
  will be scaled to have physical size \code{maxsize}.
  If \code{meansize} is given, then the
  average absolute mark value, \code{mmean=mean(abs(marks(x)))},
  will be scaled to have physical size \code{meansize}.

  The user can set the default values of these plotting parameters
  using \code{\link{spatstat.options}("par.points")}.
  
  To zoom in (to view only a subset of the point pattern at higher
  magnification), use the graphical arguments
  \code{xlim} and \code{ylim} to specify the rectangular field of view.

  The value returned by this plot function is an object of
  class \code{"symbolmap"} representing the mapping from mark values
  to graphical symbols. See \code{\link{symbolmap}}.
  It can be used to make a suitable legend,
  or to ensure that two plots use the same graphics map.
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
   plot(cells)

   plot(cells, pch=16)

   # make the plotting symbols larger (for publication at reduced scale)
   plot(cells, cex=2)

   # set it in spatstat.options
   oldopt <- spatstat.options(par.points=list(cex=2))
   plot(cells)
   spatstat.options(oldopt)

   # multitype 
   plot(lansing)

   # marked by a real number
   plot(longleaf)

   # just plot the points
   plot(longleaf, use.marks=FALSE)
   plot(unmark(longleaf)) # equivalent

   # point pattern with multiple marks
   plot(finpines)
   plot(finpines, which.marks="height")

   # controlling COLOURS of points
   plot(cells, cols="blue")
   plot(lansing, cols=c("black", "yellow", "green", 
                        "blue","red","pink"))
   plot(longleaf, fg="blue")

   # make window purple
   plot(lansing, border="purple")
   # make everything purple
   plot(lansing, border="purple", cols="purple", col.main="purple",
                 leg.args=list(col.axis="purple"))
 
   # controlling PLOT CHARACTERS
   plot(lansing, chars = 11:16)
   plot(lansing, chars = c("o","h","m",".","o","o"))

   # controlling MARK SCALE
   plot(longleaf, markscale=0.1)
   plot(longleaf, maxsize=5)
   plot(longleaf, meansize=2)

   # draw circles of diameter equal to nearest neighbour distance
   plot(cells \%mark\% nndist(cells), markscale=1, legend=FALSE)

   # inspecting the symbol map
   v <- plot(amacrine)
   v

   ## variable colours ('cols' not 'col')
   plot(longleaf, cols=function(x) ifelse(x < 30, "red", "black"))

   ## re-using the same mark scale
   a <- plot(longleaf)
   juveniles <- longleaf[marks(longleaf) < 30]
   plot(juveniles, symap=a)

   ## numerical marks mapped to symbols of fixed size with variable colour
   ra <- range(marks(longleaf))
   colmap <- colourmap(terrain.colors(20), range=ra)
   ## filled plot characters are the codes 21-25
   ## fill colour is indicated by 'bg'
   sy <- symbolmap(pch=21, bg=colmap, range=ra)
   plot(longleaf, symap=sy)

}

\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{hplot}
