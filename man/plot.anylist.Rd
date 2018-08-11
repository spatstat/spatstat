\name{plot.anylist}
\alias{plot.anylist}
\title{Plot a List of Things}
\description{
  Plots a list of things
}
\usage{
  \method{plot}{anylist}(x, \dots, main, arrange=TRUE,
   nrows=NULL, ncols=NULL, main.panel=NULL,
   mar.panel=c(2,1,1,2), hsep=0, vsep=0,
   panel.begin=NULL, panel.end=NULL, panel.args=NULL,
   panel.begin.args=NULL, panel.end.args=NULL, panel.vpad=0.2,
   plotcommand="plot",
   adorn.left=NULL, adorn.right=NULL, adorn.top=NULL, adorn.bottom=NULL,
   adorn.size=0.2, equal.scales=FALSE, halign=FALSE, valign=FALSE)
}
\arguments{
  \item{x}{
    An object of the class \code{"anylist"}.
    Essentially a list of objects.
  }
  \item{\dots}{
    Arguments passed to \code{\link{plot}} when generating each
    plot panel.
  }
  \item{main}{
    Overall heading for the plot.
  }
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
  \item{main.panel}{
    Optional. A character string, or a vector of character strings,
    giving the headings for each of the objects.
  }
  \item{mar.panel}{
    Size of the margins outside each plot panel.
    A numeric vector of length 4 giving the bottom, left, top,
    and right margins in that order. (Alternatively the vector
    may have length 1 or 2 and will be replicated to length 4).
    See the section on \emph{Spacing between plots}.
  }
  \item{hsep,vsep}{
    Additional horizontal and vertical separation between plot panels,
    expressed in the same units as \code{mar.panel}. 
  }
  \item{panel.begin,panel.end}{
    Optional. Functions that will be executed before and after each panel is
    plotted. See Details.
  }
  \item{panel.args}{
    Optional. Function that determines different plot arguments
    for different panels. See Details.
  }
  \item{panel.begin.args}{
    Optional. List of additional arguments for \code{panel.begin}
    when it is a function.
  }
  \item{panel.end.args}{
    Optional. List of additional arguments for \code{panel.end}
    when it is a function.
  }
  \item{panel.vpad}{
    Amount of extra vertical space that should be allowed for the
    title of each panel, if a title will be displayed.
    Expressed as a fraction of the height of the panel.
    Applies only when \code{equal.scales=FALSE} (the default) and
    requires that the height of each panel can be determined.
  }
  \item{plotcommand}{
    Optional.
    Character string containing the name of the command that should be
    executed to plot each panel. 
  }
  \item{adorn.left,adorn.right,adorn.top,adorn.bottom}{
    Optional. Functions (with no arguments) that will be executed to
    generate additional plots at the margins (left, right, top and/or
    bottom, respectively) of the array of plots.
  }
  \item{adorn.size}{
    Relative width (as a fraction of the other panels' widths)
    of the margin plots.
  }
  \item{equal.scales}{
    Logical value indicating whether the components
    should be plotted at (approximately) the same physical scale. 
  }
  \item{halign,valign}{
    Logical values indicating whether panels in a column
    should be aligned to the same \eqn{x} coordinate system
    (\code{halign=TRUE}) and whether panels in a row should
    be aligned to the same \eqn{y} coordinate system (\code{valign=TRUE}).
    These are applicable only if \code{equal.scales=TRUE}.
  }
}
\value{
  Null.
}
\details{
  This is the \code{plot} method for the class \code{"anylist"}.

  An object of class \code{"anylist"} represents
  a list of objects intended to be treated in the same way.
  This is the method for \code{plot}.

  In the \pkg{spatstat} package, various functions produce
  an object of class \code{"anylist"}, essentially a list of
  objects of the same kind. 
  These objects can be plotted in a nice arrangement
  using \code{plot.anylist}. See the Examples.

  The argument \code{panel.args} determines extra graphics parameters
  for each panel. It should be a function that will be called
  as \code{panel.args(i)} where \code{i} is the panel number.
  Its return value should be a list of graphics parameters that can
  be passed to the relevant \code{plot} method. These parameters
  override any parameters specified in the \code{\dots} arguments.

  The arguments \code{panel.begin} and \code{panel.end}
  determine graphics that will be plotted before and after
  each panel is plotted. They may be objects
  of some class that can be plotted
  with the generic \code{plot} command. Alternatively they
  may be functions that will be
  called as \code{panel.begin(i, y, main=main.panel[i])}
  and \code{panel.end(i, y, add=TRUE)} where \code{i} is the panel
  number and \code{y = x[[i]]}.

  If all entries of \code{x} are pixel images,
  the function \code{\link{image.listof}} is called to control
  the plotting. The arguments \code{equal.ribbon} and \code{col}
  can be used to determine the colour map or maps applied.

  If \code{equal.scales=FALSE} (the default), then the 
  plot panels will have equal height on the plot device
  (unless there is only one column of panels, in which case
  they will have equal width on the plot device). This means that the
  objects are plotted at different physical scales, by default.

  If \code{equal.scales=TRUE}, then the dimensions of the
  plot panels on the plot device will be proportional
  to the spatial dimensions of the
  corresponding components of \code{x}. This means that the
  objects will be plotted at \emph{approximately} equal physical scales.
  If these objects have very different spatial sizes,
  the plot command could fail (when it tries
  to plot the smaller objects at a tiny scale), with an error
  message that the figure margins are too large.

  The objects will be plotted at \emph{exactly} equal
  physical scales, and \emph{exactly} aligned on the device,
  under the following conditions:
  \itemize{
    \item
    every component of \code{x} is a spatial object
    whose position can be shifted by \code{\link{shift}};
    \item
    \code{panel.begin} and \code{panel.end} are either
    \code{NULL} or they are spatial objects 
    whose position can be shifted by \code{\link{shift}};
    \item
    \code{adorn.left}, 
    \code{adorn.right}, 
    \code{adorn.top} and 
    \code{adorn.bottom} are all \code{NULL}.
  }

  Another special case is when every component of \code{x} is an
  object of class \code{"fv"} representing a function.
  If \code{equal.scales=TRUE} then all these functions will 
  be plotted with the same axis scales
  (i.e. with the same \code{xlim} and the same \code{ylim}).
}
\section{Spacing between plots}{
  The spacing between individual plots is controlled by the parameters
  \code{mar.panel}, \code{hsep} and \code{vsep}.

  If \code{equal.scales=FALSE}, the plot panels are
  logically separate plots. The margins for each panel are
  determined by the argument \code{mar.panel} which becomes 
  the graphics parameter \code{mar}
  described in the help file for \code{\link{par}}.
  One unit of \code{mar} corresponds to one line of text in the margin.
  If \code{hsep} or \code{vsep} are present, \code{mar.panel}
  is augmented by \code{c(vsep, hsep, vsep, hsep)/2}.
  
  If \code{equal.scales=TRUE}, all the plot panels are drawn
  in the same coordinate system which represents a physical scale.
  The unit of measurement for \code{mar.panel[1,3]}
  is one-sixth of the greatest height of any object plotted in the same row
  of panels, and the unit for \code{mar.panel[2,4]} is one-sixth of the
  greatest width of any object plotted in the same column of panels.
  If \code{hsep} or \code{vsep} are present,
  they are interpreted in the same units as \code{mar.panel[2]}
  and \code{mar.panel[1]} respectively.
}
\seealso{
  \code{\link{contour.listof}},
  \code{\link{image.listof}},
  \code{\link{density.splitppp}}
}
\section{Error messages}{
  If the error message \sQuote{Figure margins too large}
  occurs, this generally means that one of the
  objects had a much smaller physical scale than the others.
  Ensure that \code{equal.scales=FALSE}
  and increase the values of \code{mar.panel}.
}
\examples{
 trichotomy <- list(regular=cells,
                    random=japanesepines,
                    clustered=redwood)
 K <- lapply(trichotomy, Kest)
 K <- as.anylist(K)
 plot(K, main="")

# list of 3D point patterns
 ape1 <- osteo[osteo$shortid==4, "pts", drop=TRUE]
 class(ape1)
 plot(ape1, main.panel="", mar.panel=0.1, hsep=0.7, vsep=1,
      cex=1.5, pch=21, bg='white')
}
\author{\adrian
  
  
  \rolf
  
  and \ege
  
}
\keyword{spatial}
\keyword{hplot}
