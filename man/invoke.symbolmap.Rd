\name{invoke.symbolmap}
\alias{invoke.symbolmap}
\title{
  Plot Data Using Graphics Symbol Map
}
\description{
  Apply a graphics symbol map to a vector of data values and
  plot the resulting symbols.
}
\usage{
invoke.symbolmap(map, values, x, y = NULL, \dots, add = FALSE,
               do.plot = TRUE, started = add && do.plot)
}
\arguments{
  \item{map}{
    Graphics symbol map (object of class \code{"symbolmap"}).
  }
  \item{values}{
    Vector of data that can be mapped by the symbol map.
  }
  \item{x,y}{
    Coordinate vectors for the spatial locations of the
    symbols to be plotted.
  }
  \item{\dots}{
    Additional graphics parameters.
  }
  \item{add}{
    Logical value indicating whether to add the symbols to
    an existing plot (\code{add=TRUE}) or to initialise a new
    plot (\code{add=FALSE}, the default).
  }
  \item{do.plot}{
    Logical value indicating whether to actually perform the plotting.
  }
  \item{started}{
    Logical value indicating whether the plot has already been initialised.
  }
}
\details{
  A symbol map is an association between data values and graphical symbols. 

  This command applies the symbol map \code{map} to the data 
  \code{values} and plots the resulting symbols at the locations
  given by \code{\link{xy.coords}(x,y)}.
}
\value{
  (Invisibly) the maximum diameter of the symbols, in user coordinate units.
}
\author{\adrian
  
  ,
  \rolf
  
  and \ege.
}
\seealso{
  \code{\link{plot.symbolmap}} to plot the graphics map itself.
  
  \code{\link{symbolmap}} to create a graphics map.
}
\examples{
  g <- symbolmap(range=c(-1,1),
                   shape=function(x) ifelse(x > 0, "circles", "squares"),
                   size=function(x) sqrt(ifelse(x > 0, x/pi, -x))/15,
                   bg=function(x) ifelse(x > 0, "green", "red"))
  plot(square(1), main="")
  a <- invoke.symbolmap(g, runif(10, -1, 1), runifpoint(10), add=TRUE)
  a 
}
\keyword{spatial}
\keyword{hplot}
