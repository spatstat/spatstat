\name{plot.symbolmap}
\alias{plot.symbolmap}
\title{
  Plot a Graphics Symbol Map
}
\description{
  Plot a representation of a graphics symbol map,
  similar to a plot legend.
}
\usage{
\method{plot}{symbolmap}(x, \dots, main, xlim = NULL, ylim = NULL,
                           vertical = FALSE,
                           side = c("bottom", "left", "top", "right"),
                           annotate = TRUE, labelmap = NULL, add = FALSE)
} 
\arguments{
  \item{x}{
    Graphics symbol map (object of class \code{"symbolmap"}).
  }
  \item{\dots}{
    Additional graphics arguments passed to
    \code{\link{points}}, \code{\link{symbols}} or \code{\link{axis}}.
  }
  \item{main}{
    Main title for the plot. A character string.
  }
  \item{xlim,ylim}{
    Coordinate limits for the plot.
    Numeric vectors of length 2.
  }
  \item{vertical}{
    Logical. Whether to plot the symbol map in a vertical orientation.
  }
  \item{side}{
    Character string specifying the position of the text
    that annotates the symbols.
  }
  \item{annotate}{
    Logical. Whether to annotate the symbols with labels.
  }
  \item{labelmap}{
    Transformation of the labels.
    A function or a scale factor which will be applied to
    the data values corresponding to the plotted symbols.
  }
  \item{add}{
    Logical value indicating whether to add the plot to the
    current plot (\code{add=TRUE}) or to initialise a new plot.
  }
}
\details{
  A graphics symbol map is
  an association between data values and graphical symbols. 

  This command plots the graphics symbol map itself,
  in the style of a plot legend.
}
\value{
  None.
}
\author{\adrian
  
  ,
  \rolf
  
  and \ege.
}
\seealso{
  \code{\link{symbolmap}} to create a symbol map.
  
  \code{\link{invoke.symbolmap}} to apply the symbol map to some data
  and plot the resulting symbols.
}
\examples{
  g <- symbolmap(inputs=letters[1:10], pch=11:20)
  plot(g)

  g2 <- symbolmap(range=c(-1,1),
                    shape=function(x) ifelse(x > 0, "circles", "squares"),
                    size=function(x) sqrt(ifelse(x > 0, x/pi, -x)),
                    bg = function(x) ifelse(abs(x) < 1, "red", "black"))
  plot(g2, vertical=TRUE, side="left", col.axis="blue", cex.axis=2)
}
\keyword{spatial}
\keyword{hplot}
