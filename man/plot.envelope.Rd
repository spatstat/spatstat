\name{plot.envelope}
\alias{plot.envelope}
\title{Plot a Simulation Envelope}
\description{
  Plot method for the class \code{"envelope"}.
}
\usage{
 \method{plot}{envelope}(x, \dots, shade=c("hi", "lo"))
}
\arguments{
  \item{x}{
    An object of class \code{"envelope"}, containing the variables to be plotted
    or variables from which the plotting coordinates can be computed.
  }
  \item{\dots}{
    Extra arguments passed to \code{\link{plot.fv}}.
  }
  \item{shade}{
    An index that identifies two columns of \code{x}. When the corresponding
    curves are plotted, the region between the curves will be 
    shaded in light grey. 
  }
}
\value{
  Either \code{NULL}, or a data frame giving the meaning of the
  different line types and colours.
}
\details{
  This is the \code{plot} method for the class \code{"envelope"}
  of simulation envelopes. Objects of this class are
  created by the command \code{\link{envelope}}.

  This plot method is identical to \code{\link{plot.fv}} except that
  its default behaviour is to shade the region
  between the upper and lower envelopes in a light grey colour.
  
  To suppress the shading and plot the upper and lower envelopes
  as curves, set \code{shade=NULL}.

  To change the colour of the shading, use the argument \code{shadecol}
  which is passed to \code{\link{plot.fv}}. 
  
  See \code{\link{plot.fv}} for further information on how to
  control the plot.
}
\examples{
   data(cells)
   E <- envelope(cells, Kest, nsim=19)
   plot(E)
   plot(E, sqrt(./pi) ~ r)
}
\seealso{
  \code{\link{envelope}},
  \code{\link{plot.fv}}
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{hplot}
 
