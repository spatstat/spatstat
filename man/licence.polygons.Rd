\name{licence.polygons}
\alias{licence.polygons}
\alias{license.polygons}
\title{Restricted Licence Conditions for Polygon Calculations}
\description{
  Prints a summary of the restricted licence conditions
  for polygon calculations using the package \pkg{gpclib}.
}
\usage{
  licence.polygons()
  license.polygons()
}
\value{
  Null.
}
\details{
  \emph{The functions named above simply print the following message.}
  
Polygon geometry computations in \pkg{spatstat}
(such as computing the intersection or union of two polygons)
currently depend on the package \pkg{gpclib}.

'gpclib' has a restricted licence which explicitly forbids commercial use.

Accordingly, in \pkg{spatstat} (versions \code{1.17-1} and later) 
the code which depends on \pkg{gpclib} can be either enabled or disabled.

By default, \pkg{spatstat} does NOT use \pkg{gpclib}.
In this mode, polygon geometry computations will be executed by
converting the polygons to binary pixel images (`masks'), so that 
the result of the intersection or union will be a binary pixel mask.

If you are covered by the \pkg{gpclib} licence, you can enable the use
of \pkg{gpclib} inside \pkg{spatstat} by typing
\code{spatstat.options(gpclib=TRUE)}.

To check whether \pkg{gpclib} is currently enabled or disabled, type 
\code{spatstat.options("gpclib")}.
}
\seealso{
  \code{\link{spatstat.options}}
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
